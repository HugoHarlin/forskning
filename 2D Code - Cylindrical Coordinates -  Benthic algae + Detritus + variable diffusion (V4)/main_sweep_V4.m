
%% Implementation details and a clean slate
close all
clear
clc
tic

% Version 4 - detritus + variable diffusion coefficients.

% This code implements a system of equations based on the work by Jäger et
% al. (2010) with an additional benthic algal layer at the bottom,
% in two spatial dimensions, using a finite volume method and Gauss divergence
% theorem to solve the system of equations.
% ode15s is used to iterate in time.
%
% Plankton nutrient density stochiometry q is constant, and the lake topography
% is described by an upsidedown cone where the bottom edge is sloped so
% that the lake is deepest in the middle and shallow at the shore. The slope of the
% bottom is governed by an exponent alpha, and half the lake is modelled so that the
% deepest edge corresponds the middle of the lake.
%
% The matrices in this script are indexed so that origo is centered at the center of the lake surface,
% with the z-axis facing downwards and the x-axis facing the lake shore.
%
% Written by Hugo Harlin 2019-2020

%% loop over different background turbidies
%parfor p = 1:3 % looping over diffusion coeff

%% Loop over different turbidities
for p = 1 % :3 %2:3 %1:2
    
    for depth_var = 1% :3 % looping over thermocline depth
        % for x_res_index = [0,1]
        %for therm_depth = [2,5,10,15]
        kbg_index = p;
        
        tic
        
        %% all parameters are organized in a struct pp for implementation convenience
        pp = struct;
        
        %% Lake topology and Mesh
        % Quantities relating to system size
        pp.Xn   = 15;  % Number of grid-points (width)
        pp.Zn   = 15;  % Number of grid-points (depth)
        pp.Lmin = 0.0001; % Minimum lake depth (depth at land-water interface) [m]
        pp.Lmax = 20;  % Maximum lake depth [m]
        pp.W    = 20;  % Lake radius [m]
        alpha_vec = [1, 1.5]; % slopes being looped over.
        alpha_index = 1;
        pp.alpha = alpha_vec(alpha_index);  % Exponent governing the slope of the lake bottom
        
        pp.stratified = 0; % if true, the lake diffusion coefficients are set so that the diffusion coefficient reflect a thermally stratified
        if(pp.stratified)
            pp.increased_x_res = 0;
            pp.increased_z_res = 0;
            % lake with an epilimnion at 5 meters, (3 meters wide transition)
            depth_vec = [2, 5, 10, 15];
            
            pp.thermocline_depth = depth_vec(depth_var); % vertical depth at which the water transitions from a well mixed state, the upper boundary of the thermocline [m].
            pp.thermocline_thickness = 3; % the thickness of the thermocline layer [m]
            pp.diff_above_thermocline = 1000; % [day^-1]
            pp.diff_in_thermocline = 1; % [day^-1]
            pp.diff_below_thermocline = 10;   % [day^-1]
        end
        
        % Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
        % at the center of the lake (slope = alpha* (Lmin - Lmax)/W).
        % (0,0) is placed at the center of the lake at the surface, y-dim is facing
        % downward and x-dim is facing towards the lake edge
        
        pp.X = zeros(pp.Zn, pp.Xn); % Mesh-spacing in x-dimension
        pp.Z = zeros(pp.Zn, pp.Xn); % Mesh-spacing in y-dimension
        
        for ii=1:pp.Zn
            pp.X(ii,:) = [0:pp.W/(pp.Xn-1):pp.W]; % even spacing of the grid in x-dimension
        end
        
        % the grid is compressed in y-dimension, with depth Lmin at the
        % shore and Lmax at the center of the lake.
        for jj = 1:pp.Xn
            for iii=1:pp.Zn
                pp.Z(iii,jj) = (iii-1)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(iii,jj)/pp.W).^(pp.alpha));
            end
        end
        
        % Here the grid is refined in the z-direction to increase resolution at
        % the thermocline.
        counter = 0; % counter is used when appending the new rows in the grid
        
        if(pp.stratified)
            
            if(pp.increased_z_res) % horizontal grid refinement (more rows) where the thermocline meets the bottom
                jjjj = 1:pp.Xn;
                for iiii = pp.thermocline_depth+1:  pp.thermocline_depth +  pp.thermocline_thickness +1
                    new_row_1 = (iiii-1+0.25)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(iiii,jjjj)/pp.W).^(pp.alpha));
                    new_row_2 = (iiii-1+0.50)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(iiii,jjjj)/pp.W).^(pp.alpha));
                    new_row_3 = (iiii-1+0.75)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(iiii,jjjj)/pp.W).^(pp.alpha));
                    
                    pp.Z = [ pp.Z(1:iiii + counter,jjjj); new_row_1; new_row_2; new_row_3; pp.Z(((iiii+ counter+1):end),jjjj)];
                    counter = counter + 3; % three new rows are added.
                end
                
                % adding rows to pp.X so that the dimensions don't differ.
                pp.X = [pp.X; pp.X(1:3* (pp.thermocline_thickness+1),:)];
                pp.Zn = pp.Zn + 3* (pp.thermocline_thickness+1); % updating the count of number of rows in our grid p.Zn
            end
            
            if(pp.increased_x_res) % horizontal grid refinement (more columns) where the thermocline meets the bottom.
                x_strat_top = pp.W*((pp.thermocline_depth - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the top of the thermocline intersects with the bottom.
                x_strat_bottom = pp.W*((pp.thermocline_depth + pp.thermocline_thickness - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the bottom of the thermocline intersects with the bottom.
                x_strat_top = ceil(x_strat_top); % by rounding up to the nearest whole number, ensuring that the variable is compatible as an index and that we dont miss any grid refinement
                x_strat_bottom = floor(x_strat_bottom); % where the thermocline is.
                
                
                for ii=1:1:pp.Zn
                    pp.X(ii,:) = [0:pp.W/(pp.Xn-1):pp.W]; % even spacing of the grid in x-dimension
                end
                
                % the grid is compressed in y-dimension, with depth Lmin at the
                % shore and Lmax at the center of the lake.
                tt=1:pp.Zn;
                
                temp_index = x_strat_bottom +1 ; % index used for concatenating existing Z-matrix with new rows
                
                for j = x_strat_bottom:x_strat_top-1
                    %for j = [1]
                    num_new_columns = 1;
                    new_columns = zeros(pp.Zn,num_new_columns); % new columnds for pp.Z matrix
                    temp_X = zeros(pp.Zn,num_new_columns); % new colums for pp.X matrix
                    
                    for k = 1:num_new_columns
                        new_columns(:,k) = (tt-1)'.*(pp.Lmax/(pp.Zn-1)).*(1 + (pp.Lmin/pp.Lmax -1).*((pp.X(tt, temp_index)+k/( num_new_columns +1) )/pp.W).^(pp.alpha));
                        temp_X(:,k) = pp.X(tt,temp_index) + k/(num_new_columns+1);
                    end
                    
                    pp.X  = [pp.X(:,1:temp_index)  temp_X  pp.X(:,temp_index+1:end)];
                    
                    pp.Z = [pp.Z(:,1:temp_index) new_columns pp.Z(:,temp_index+1:end)]; % adding in the new columns in the Z-matrix
                    temp_index = temp_index +num_new_columns+1; % adjusting the index where the new rows are to be added.s
                end
                
                pp.Xn = pp.Xn + num_new_columns* (j-x_strat_bottom+1);
            end
            
        end
        
        
        % coordinates of the center of each mesh quadrilateral
        X_vol = zeros(pp.Zn-1, pp.Xn-1);
        Z_vol = zeros(pp.Zn-1, pp.Xn-1);
        
        for j=1:pp.Xn-1
            for ttt =1:pp.Zn-1
                X_vol(ttt,j) = (pp.X(ttt,j) + pp.X(ttt+1,j) + pp.X(ttt,j+1) + pp.X(ttt+1,j+1) ) /4;
                Z_vol(ttt,j) = (pp.Z(ttt,j) + pp.Z(ttt+1,j) + pp.Z(ttt,j+1) + pp.Z(ttt+1,j+1) ) /4;
            end
        end
        pp.X_vol = X_vol;
        pp.Z_vol = Z_vol;
        
        pp.volumes_cyl = vol_areas_cyl_3d_fn(pp); % volumes of each grid cell in the bulk of the lake
        pp.Area_bottom_cyl = Area_bottom_cyl_fn(pp); % area of each grid element on the bottom
        
        %% Model Parameters
        
        %%% Diffusion coefficients are defined below %%%
        pp.I0 = 300;    % Light intensity at the surface [micro-mol photons m^-2 s^-1]
        pp.kA = 0.0003; % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
        pp.kD = 0.0003; % Specific light-attenuation coefficient of detritus [m^2 mg C^-1]
        pp.kB =  0.00003; % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
        %kgb_vec = [0 0.2 0.4 0.8 2.0];
        kgb_vec = [0 0.2 0.8 2.0];
         kbg_index = 2;
        pp.kbg =  kgb_vec(kbg_index);  % Background light-attenuation coefficient [m^-1]
        pp.lbg_A = 0.08; % Specific algal maintenance respiration losses [day^-1]
        pp.Ad = 0.02;    % algal death rate [day^-1]
        pp.M = 1.5;      % Half-saturation constant of algal nutrient uptake [mg P m^-3]
        pp.M_benth = 1.5; % Half-saturation constant of benthic algae nutrient uptake [mgPm^-3]
        pp.Gmax = 1.08;   % Maximum specific phytoplankton production rate [day^-1]
        pp.Gmax_benth = 1.08; % Maximum specific benthic algae production rate [day^-1]
        pp.H = 120.0;         % Half-saturation constant of light-dependent algal production [micro-mol photons m^-2 s^-1]
        pp.H_benth = 120.0;   % Half-saturation constant of light-dependent benthic algal production [micro-mol photons m^-2 s^-1]
        pp.q = 0.0244;        % Algal nutrient quota, Redfield ratio [mgP/mgC]
        pp.q_benth =  0.0244; % Benthic algae nutrient quota, Redfield ratio [mgP/mgC]
        pp.lbg_benth = 0.1;   % Specific benthic algae maintenance respiration losses [day^-1]
        remin_vec = [0, 0.01,0.02, 0.05, 0.1]; % 0.02 is used as default value
        remin_index = 2;
        pp.r = remin_vec(remin_index);  % Specific mineralization rate of sedimented nutrients [day^-1]
        pp.vA = 0.1;          % Algal sinking speed [m day^-1]
        pp.vD = 0.25;         % detritus sinking speed [m day^-1]
        pp.Dbg = 0.02;        % remineralization of detritus in the water column [day^-1]
        pp.resus_Max = 0.02;  % Maximum resuspension rate of the detritus in the sediment [day^-1]
        pp.resus_H = 8;       % half saturation constant of the resuspension rate functional response [m^2 day^-1]
        pp.death_rate = 1;    % coefficient governing the proportion of sinking algae at the bottom that dies.
        % 0 = no death. 1 = all algae that would have sunk through the sediment dies.
        pp.benth_recycling = 0.5; % range: [0,1]. Governs the portion of respired nutrients that are released as dissolved nutrients.
        % the rest is bound in particulate matter in the sediment.
        
        pp.constant_resuspension = 1;
        
        if(pp.constant_resuspension)
            resus_vec = [0, 0.1, 0.2]; % 0.2 is used as default value
            resus_index = 2;
            pp.resus = resus_vec(resus_index);  % resuspension rate of the detritus from the sediment. [day^-1]
            
        else
            pp.resus_depth = ones(1,pp.Xn-1); % Resuspension rate that varies with depth
            max_resus = 0.2;
            min_resus = 0.002;
            i=(0:pp.Xn-2);
            pp.resus_depth(i+1) = min_resus +  i*(max_resus- min_resus)/(pp.Xn-2);
        end
        %%% periodic mixing
        pp.periodic_mixing = 0; %set to true if periodic mixing is desired
        pp.mixing_frequency = 365; % time between each mixing event [days]
        pp.max_sim_time = 1e+4; %10e+6; % upper bound of the total simulation time if steady state is not reached before then. [days]
        
        %% looping over the diffusion coefficients
        % for diff_index = 1:3 %,1000]
        %for diff_max_z = [1,10,100] %,1000]
        
        %% Diffusion Coefficients
        dx = zeros(pp.Zn-1,pp.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
        dz = zeros(pp.Zn-1,pp.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
        
        if(~pp.stratified)
            diff_vec =[0,1,10,100];
            diff_max_x =diff_vec(2);
           % diff_max_x = diff_vec(diff_index); % horizontal diffusion coefficient at the surface
           % diff_max_x = 0;
            %diff_max_z = 10; % vertical diffusion coefficient at the surface
            diff_min_x = diff_max_x; % horizontal diffusion coefficient at the bottom
            diff_min_z = diff_max_x; % vertical diffusion coefficient at the bottom
            diff_max_z = diff_max_x;
            %diff_min_z = diff_max_z; % vertical diffusion coefficient at the bottom
            
            
            % Linearly decreasing diffusion coefficients with the lake depth.
            for j = 1:pp.Xn-1
                for ty = 1:pp.Zn-1
                    dx(ty,j) =  diff_min_x + (diff_max_x-diff_min_x)*(1- pp.Z_vol(ty,j)/pp.Lmax);
                    dz(ty,j) =  diff_min_z + (diff_max_z-diff_min_z)*(1- pp.Z_vol(ty,j)/pp.Lmax);
                end
            end
            
            pp.dx = dx;
            pp.dz = dz;
        end
        
        %pp.dx = diff_max_z.*dx;
        %pp.dz = diff_max_z.*dz;
        
        if(pp.stratified)
            %%%%% manual setting of the diffusion coefficients %%%%
            for x = 1:pp.Xn-1
                for z = 1:pp.Zn-1
                    if(pp.Z_vol(z,x) < pp.thermocline_depth)
                        dx(z,x) =  pp.diff_above_thermocline;
                        dz(z,x) =  pp.diff_above_thermocline;
                    elseif(pp.thermocline_depth < pp.Z_vol(z,x) &&  pp.Z(z,x) < pp.thermocline_depth + pp.thermocline_thickness)
                        dx(z,x) = pp.diff_in_thermocline;
                        dz(z,x) = pp.diff_in_thermocline;
                    else
                        dx(z,x) =  pp.diff_below_thermocline;
                        dz(z,x) =  pp.diff_below_thermocline;
                    end
                end
            end
        end
        
        
        pp.dx = dx;
        pp.dz = dx;
        
        %% variable resuspension coefficient
        % A type II functional response function is used to
        %p.resus_depth
        
        %% Pre-evaluation of transform derivatives and creation of stiffness matrix
        [dX_dXi_preCalc, dX_dEta_preCalc, dY_dXi_preCalc, dY_dEta_preCalc] = preCalc_Derivatives(pp);
        pp.dX_dXi_preCalc   = dX_dXi_preCalc;
        pp.dX_dEta_preCalc  = dX_dEta_preCalc;
        pp.dZ_dXi_preCalc   = dY_dXi_preCalc;
        pp.dZ_dEta_preCalc  = dY_dEta_preCalc;
        
        %% Stiffness matrix & light integration matrix
        % the stiffness matrix is the matrix M in the linearized system y' = My,
        % where y is a column vector.
        pp.S = Stiffness_matrix(pp);
        pp.I_matrix = I_matrix_fn(pp);
        
        %% Inital Conditions
        A0  = 1.0*ones(pp.Zn-1, pp.Xn-1);    % Initial Algal carbon density [mg C m^-3]
        Rd0 = 10*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of dissolved nutrients [mg P m^-3]
        D0  = 1.000*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of detritus [mg P m^-3]
        Rs0 = 1.0*ones(1, pp.Xn-1);        % initial concentration of sediment nutrient density [mg P m^-2]
        B0  = 1.00*ones(1, pp.Xn-1);        % initial concentration of benthic algal density [mg C m^-2]
        
       % A0(end,:) = 1;
        %% calculation of total nutrient content at t=0, (to be conserved)
        n_algae_0 = pp.volumes_cyl.*A0*pp.q;
        n_dissolved_0 = pp.volumes_cyl.*Rd0;
        n_detritus_0 = pp.volumes_cyl.*D0;
        n_sediment_0 = Rs0.*pp.Area_bottom_cyl;
        n_benthic_0 = B0.*pp.q_benth.*pp.Area_bottom_cyl;
        
        pp.ntot_algae_0  = sum(sum(n_algae_0));
        pp.ntot_dissolved_0 = sum(sum(n_dissolved_0));
        pp.ntot_detritus_0 = sum(sum(n_detritus_0));
        pp.ntot_sed_0 = sum(n_sediment_0);
        pp.ntot_benthic_0 = sum(n_benthic_0);
        
        ntot_0 =  pp.ntot_algae_0 + pp.ntot_dissolved_0 + pp.ntot_detritus_0 + pp.ntot_sed_0 + pp.ntot_benthic_0;
        pp. ntot_0 =  ntot_0;
        
        %% state variables
        A = A0'; % transpose of matrices in order to use colon notation to reshape to vector form.
        Rd = Rd0';
        D = D0';
        Rs = Rs0';
        B = B0';
        y0 = [A(:); Rd(:); D(:); Rs(:); B(:)];
        
        %% Simulation of model
        
        
        ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events', @(t,y) eventfun_V4(t,y,pp), 'NonNegative',(1:length(y0)));
        
        if(~pp.periodic_mixing) % no periodic mixing
            %tend = 2000; % end simulation time
            %t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
            %t_span = [1:tend];
            %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events',@(t,y) eventfun_V4(t,y,p) , 'NonNegative',(1:length(y0)));
            %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'NonNegative',find(y0));
            [t,Y_t] = ode15s( @(t,Y) rhs_function_V4(t,Y,pp), [0, inf] , y0 , ode_opts);
            %[t,Y_t] = ode15s( @(t,Y) dAdt_efficient_correct_V4_detritus(t,Y,p) , t_span , y0 , ode_opts);
            
        else % periodic mixing
            t_var = 0; % variable keeping track of the total lapsed simulation time
            res_total = []; % state variables for all time states of the total simulation, analog of Y_t for a single ode15s run
            time_total = []; % this vector will contain the simulation time for the total simulation, analog to the t vector if ode15s is only run once.
            
            while(t_var <  pp.max_sim_time)
                [t,Y_t] = ode15s( @(t,Y) rhs_function_V4(t,Y,pp), [ t_var,  t_var+pp.mixing_frequency] , y0 , ode_opts);
                
                % saving timesteps and state variables
                time_total = [time_total;  t];
                res_total = [res_total; Y_t];
                
                % updating the timekeeping variable
                t_var = t_var +pp.mixing_frequency;
                
                [xx,isterm,dir] = eventfun_V4(t(end),Y_t(end,:)',pp);
                
                if(xx >0) % steady state is not reached and we have arrived at the mixing event
                    %%%%%%%%%%% mixing event %%%%%%%%%%%%
                    % Extracting results
                    A  = Y_t(end,1:(pp.Xn-1)*(pp.Zn-1));
                    Rd = Y_t(end,(pp.Xn-1)*(pp.Zn-1)+1 : 2*(pp.Xn-1)*(pp.Zn-1));
                    D  = Y_t(end,2*(pp.Xn-1)*(pp.Zn-1)+1 : 3*(pp.Xn-1)*(pp.Zn-1));
                    Rs = Y_t(end,3*(pp.Xn-1)*(pp.Zn-1) +1 : 3*(pp.Xn-1)*(pp.Zn-1) + (pp.Xn-1));
                    B  = Y_t(end,3*(pp.Xn-1)*(pp.Zn-1) + (pp.Xn-1)+1: end);
                    
                    % Reshaping the vectors into matrices
                    A  = reshape(A, [pp.Xn-1, pp.Zn-1]);
                    Rd = reshape(Rd, [pp.Xn-1, pp.Zn-1]);
                    D  = reshape(D, [pp.Xn-1, pp.Zn-1]);
                    A  = A';
                    Rd = Rd';
                    D  = D';
                    
                    %Calculating nutrient content in each grid cell at time t
                    n_algae     = pp.volumes_cyl.*A;
                    n_dissolved = pp.volumes_cyl.*Rd;
                    n_detritus  = pp.volumes_cyl.*D;
                    
                    %calculating concentration if well mixed
                    well_mixed_A =  sum(sum(n_algae)) / sum(sum(pp.volumes_cyl));
                    well_mixed_Rd =  sum(sum(n_dissolved)) / sum(sum(pp.volumes_cyl));
                    well_mixed_D = sum(sum(n_detritus)) / sum(sum(pp.volumes_cyl));
                    
                    % new concentration when well mixed
                    A =  well_mixed_A.*ones(pp.Zn-1, pp.Xn-1);  % phytoplankton
                    Rd = well_mixed_Rd.*ones(pp.Zn-1, pp.Xn-1);  % disolved nutrients
                    D = well_mixed_D.*ones(pp.Zn-1, pp.Xn-1);   % detritus
                    
                    y0 = [A(:); Rd(:); D(:); Rs(:); B(:)]; % state variable input for the next function call of ode15s, with the bulk of the lake being well mixed.
                    
                else
                    t_var = pp.max_sim_time +1; % we exit the loop if steady state is reached.
                end
            end
            Y_t = res_total;
            t = time_total;
            
        end
          
        %% recording of simulation time & saving workspace
        pp.simTime = toc;
        alpha_str = ["1_0","1_5"];
        resus_str = ["0_0", "0_1","0_2"];
        kbg_str = ["0_2", "0_8", "2_0"];
        % kbg_str = ["0", "0_2", "0_4", "0_8", "2_0"];
        remin_str = ["0_01","0_05","0_1"];
        
        
        file_name = "2D_benthic_results_V4_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) ;
        if(pp.stratified)
            file_name = file_name +  "_Stratified_therm_depth_" + num2str(pp.thermocline_depth);
            if(pp.increased_x_res)
                file_name = file_name + "_increased_x_res";
            end
        else
            file_name = file_name  +"_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z);
        end
        
        if(pp.periodic_mixing)
            file_name = file_name +  "_periodic_mixing_" + num2str(pp.mixing_frequency) + "_end_sim_time_" + num2str( pp.max_sim_time);
        end
        
        if(0)
            %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_resuspRate_"+ resus_str(resus_index) ;
            %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_reminrate_"+ remin_str(remin_index) ;
            %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_init_dissolved_P_1000" ;
            file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_linearly_receding_resusp";
        end
        
        parsave(file_name,pp,Y_t,t,simTime,y0);
        
        
        % end
    end
end

function parsave(fname, p,Y_t,t,simTime,y0)
save(fname, 'p', 'Y_t','t','simTime','y0');
end
