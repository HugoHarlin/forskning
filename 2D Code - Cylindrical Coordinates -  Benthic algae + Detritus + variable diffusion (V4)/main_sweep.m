
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
for p = 2:5
    
    for x_res_index = [0,1]
        for therm_depth = [2,5,10,15]
            kbg_index = p;
            
            tic
            
            %% all parameters are organized in a struct pp for implementation convenience
            pp = struct;
            
            %% Lake topology and Mesh
            % Quantities relating to system size
            pp.Xn   = 21;  % Number of grid-points (width)
            pp.Zn   = 21;  % Number of grid-points (depth)
            pp.Lmin = 0.1; % Minimum lake depth (depth at land-water interface) [m]
            pp.Lmax = 20;  % Maximum lake depth [m]
            pp.W    = 20;  % Lake radius [m]
            alpha_vec = [1, 1.5]; % slopes being looped over.
            alpha_index = 2;
            pp.alpha = alpha_vec(alpha_index);  % Exponent governing the slope of the lake bottom
            
            stratified = true; % if true, the lake diffusion coefficients are set so that the diffusion coefficient reflect a thermally stratified
            increased_x_res = 1;
            increased_z_res = 0;
            % lake with an epilimnion at 5 meters, (3 meters wide transition)
            thermocline_depth = therm_depth; % vertical depth at which the water transitions from a well mixed state, the upper boundary of the thermocline [m].
            thermocline_thickness = 3; % the thickness of the thermocline layer [m]
            diff_above_thermocline = 1000; % [day^-1]
            diff_in_thermocline = 1; % [day^-1]
            diff_below_thermocline = 10;   % [day^-1]
            
            
            % Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
            % at the center of the lake (slope = alpha* (Lmin - Lmax)/W).
            % (0,0) is placed at the center of the lake at the surface, y-dim is facing
            % downward and x-dim is facing towards the lake edge
            
            pp.X = zeros(pp.Zn, pp.Xn); % Mesh-spacing in x-dimension
            pp.Z = zeros(pp.Zn, pp.Xn); % Mesh-spacing in y-dimension
            
            for i=1:1:pp.Zn
                pp.X(i,:) = [0:pp.W/(pp.Xn-1):pp.W]; % even spacing of the grid in x-dimension
            end
            
            % the grid is compressed in y-dimension, with depth Lmin at the
            % shore and Lmax at the center of the lake.
            for j = 1:pp.Xn
                for i=1:pp.Zn
                    pp.Z(i,j) = (i-1)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(i,j)/pp.W).^(pp.alpha));
                end
            end
            
            % Here the grid is refined in the z-direction to increase resolution at
            % the thermocline.
            counter = 0; % counter is used when appending the new rows in the grid
            
            if(stratified)
                
                if(increased_z_res) % horizontal grid refinement (more rows) where the thermocline meets the bottom
                    j = 1:pp.Xn;
                    for i = thermocline_depth+1:  thermocline_depth +  thermocline_thickness +1
                        new_row_1 = (i-1+0.25)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(i,j)/pp.W).^(pp.alpha));
                        new_row_2 = (i-1+0.50)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(i,j)/pp.W).^(pp.alpha));
                        new_row_3 = (i-1+0.75)*(pp.Lmax/(pp.Zn-1))*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(i,j)/pp.W).^(pp.alpha));
                        
                        test = 1;
                        pp.Z = [ pp.Z(1:i + counter,j); new_row_1; new_row_2; new_row_3; pp.Z(((i+ counter+1):end),j)];
                        counter = counter + 3; % three new rows are added.
                    end
                    
                    % adding rows to pp.X so that the dimensions don't differ.
                    pp.X = [pp.X; pp.X(1:3* (thermocline_thickness+1),:)];
                    pp.Zn = pp.Zn + 3* (thermocline_thickness+1); % updating the count of number of rows in our grid p.Zn
                end
                
                if(increased_x_res) % horizontal grid refinement (more columns) where the thermocline meets the bottom.
                    x_strat_top = pp.W*((thermocline_depth - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the top of the thermocline intersects with the bottom.
                    x_strat_bottom = pp.W*((thermocline_depth + thermocline_thickness - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the bottom of the thermocline intersects with the bottom.
                    x_strat_top = ceil(x_strat_top); % by rounding up to the nearest whole number, ensuring that the variable is compatible as an index and that we dont miss any grid refinement
                    x_strat_bottom = floor(x_strat_bottom); % where the thermocline is.
                    
                    
                    for i=1:1:pp.Zn
                        pp.X(i,:) = [0:pp.W/(pp.Xn-1):pp.W]; % even spacing of the grid in x-dimension
                    end
                    
                    % the grid is compressed in y-dimension, with depth Lmin at the
                    % shore and Lmax at the center of the lake.
                    i=1:pp.Zn;
                    
                    temp_index = x_strat_bottom +1 ; % index used for concatenating existing Z-matrix with new rows
                    
                    for j = x_strat_bottom:x_strat_top-1
                        %for j = [1]
                        num_new_columns = 1;
                        new_columns = zeros(pp.Zn,num_new_columns); % new columnds for pp.Z matrix
                        temp_X = zeros(pp.Zn,num_new_columns); % new colums for pp.X matrix
                        
                        for k = 1:num_new_columns
                            new_columns(:,k) = (i-1)'.*(pp.Lmax/(pp.Zn-1)).*(1 + (pp.Lmin/pp.Lmax -1).*((pp.X(i, temp_index)+k/( num_new_columns +1) )/pp.W).^(pp.alpha));
                            temp_X(:,k) = pp.X(i,temp_index) + k/(num_new_columns+1);
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
                for i =1:pp.Zn-1
                    X_vol(i,j) = (pp.X(i,j) + pp.X(i+1,j) + pp.X(i,j+1) + pp.X(i+1,j+1) ) /4;
                    Z_vol(i,j) = (pp.Z(i,j) + pp.Z(i+1,j) + pp.Z(i,j+1) + pp.Z(i+1,j+1) ) /4;
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
            pp.kB = 0.00003; % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
            kgb_vec = [0 0.2 0.4 0.8 2.0];
            % kbg_index = 2;
            pp.kbg = kgb_vec(kbg_index);  % Background light-attenuation coefficient [m^-1]
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
            remin_vec = [0.01,0.02, 0.05, 0.1]; % 0.02 is used as default value
            remin_index = 3;
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
            
            resus_vec = [0, 0.1, 0.2]; % 0.2 is used as default value
            resus_index = 3;
            pp.resus = resus_vec(resus_index);  % resuspension rate of the detritus from the sediment. [day^-1]
            
            pp.resus_depth = ones(1,pp.Xn-1); % Resuspension rate that varies with depth
            max_resus = 0.2;
            min_resus = 0.002;
            ii=(0:pp.Xn-2);
            pp.resus_depth(ii+1) = min_resus +  ii*(max_resus- min_resus)/(pp.Xn-2);
            
            %% looping over the diffusion coefficients
            %for diff_max_x = [1,10,100,1000]
            for diff_max_z = [1,10,100] %,1000]
                
                %% Diffusion Coefficients
                dx = zeros(pp.Zn-1,pp.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
                dz = zeros(pp.Zn-1,pp.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
                
                %     if(~stratified)
                %         diff_vec =[1,10,100,1000];
                %         diff_max_x =diff_vec(p);
                %         %diff_max_x = 100; % horizontal diffusion coefficient at the surface
                %         %diff_max_z = 100; % vertical diffusion coefficient at the surface
                %         diff_min_x = diff_max_x; % horizontal diffusion coefficient at the bottom
                %         diff_min_z = diff_max_z; % vertical diffusion coefficient at the bottom
                %
                %
                %         % Linearly decreasing diffusion coefficients with the lake depth.
                %         for j = 1:pp.Xn-1
                %             for i = 1:pp.Zn-1
                %                 dx(i,j) =  diff_min_x + (diff_max_x-diff_min_x)*(1- pp.Z_vol(i,j)/pp.Lmax);
                %                 dz(i,j) =  diff_min_z + (diff_max_z-diff_min_z)*(1- pp.Z_vol(i,j)/pp.Lmax);
                %             end
                %         end
                %
                %pp.dx = dx;
                %pp.dz = dz;
                %     end
                
                %pp.dx = diff_max_z.*dx;
                %pp.dz = diff_max_z.*dz;
                
                if(stratified)
                    %%%%% manual setting of the diffusion coefficients %%%%
                    for x = 1:pp.Xn-1
                        for z = 1:pp.Zn-1
                            if(pp.Z_vol(z,x) < thermocline_depth)
                                dx(z,x) =  diff_above_thermocline;
                                dz(z,x) =  diff_above_thermocline;
                            elseif(thermocline_depth < pp.Z_vol(z,x) &&  pp.Z(z,x) < thermocline_depth + thermocline_thickness)
                                dx(z,x) = diff_in_thermocline;
                                dz(z,x) = diff_in_thermocline;
                            else
                                dx(z,x) =  diff_below_thermocline;
                                dz(z,x) =  diff_below_thermocline;
                            end
                        end
                    end
                end
                
                
                pp.dx = dx;
                pp.dz = dx;
                
            end
            
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
            Rd0 = 10.000*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of dissolved nutrients [mg P m^-3]
            D0  = 1.000*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of detritus [mg P m^-3]
            Rs0 = 1.0*ones(1, pp.Xn-1);        % initial concentration of sediment nutrient density [mg P m^-2]
            B0  = 1.00*ones(1, pp.Xn-1);        % initial concentration of benthic algal density [mg C m^-2]
            %A0(1,:)= 1.0;
            %  A0(:,3)= 1.0;
            %  A0(:,5)= 1.0;
            %  A0(:,7)= 1.0;
            %  A0(:,9)= 1.0;
            %  A0(:,11)= 1.0;
            %  A0(:,13)= 1.0;
            %  A0(:,15)= 1.0;
            % A0(1,3)= 1.5;
            % A0(1,1)= 1;
            % A0(5,5) = 0;
            % A0(2,:) = 1;
            % A0(3,:) = 1.1;
            %A0(1,:) = 10./p.volumes_cyl(1,:);
            %A0(15,15)=10;
            %A0(:,p.Xn-4) = 10;
            % A0(1,1) = 1;
            
            % Rd0(1,1) = 1;
            % Rd0(1,2) = 1.1;
            %Rd0(2,2) = 1.1;
            % Rd0(2,2) = 2;
            %Rd0(6,7) = 1;
            %Rd0(16,16) = 5;
            
            %Rd0(6,7) = 1;
            
            %Rs0(1) = 1;
            
            %B0(1) = 0;
            
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
            %tend = 2000; % end simulation time
            %t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
            %t_span = [1:tend];
            %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events',@(t,y) eventfun_V4(t,y,p) , 'NonNegative',(1:length(y0)));
            %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'NonNegative',find(y0));
            ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events', @(t,y) eventfun_V4(t,y,pp), 'NonNegative',(1:length(y0)));
            [t,Y_t] = ode15s( @(t,Y) rhs_function_V4(t,Y,pp), [0, inf] , y0 , ode_opts);
            %[t,Y_t] = ode15s( @(t,Y) dAdt_efficient_correct_V4_detritus(t,Y,p) , t_span , y0 , ode_opts);
            
            %% recording of simulation time & saving workspace
            simTime = toc;
            alpha_str = ["1_0","1_5"];
            resus_str = ["0_0", "0_1","0_2"];
            %kbg_str = ["0_4", "0_8"];
            kbg_str = ["0", "0_2", "0_4", "0_8", "2_0"];
            remin_str = ["0_01","0_05","0_1"];
            
            
            if(stratified)
                file_name = "2D_benthic_results_V4_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_Stratified_therm_depth_" + num2str(therm_depth);
                if(increased_x_res)
                    file_name = file_name + "_increased_x_res";
                end
            else
                %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_resuspRate_"+ resus_str(resus_index) ;
                %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_reminrate_"+ remin_str(remin_index) ;
                %file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_init_dissolved_P_1000" ;
                file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_x) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_alpha_" + alpha_str(alpha_index) + "_kbg_"+ kbg_str(kbg_index) +  "_linearly_receding_resusp";
            end
            parsave(file_name,pp,Y_t,t,simTime,y0);
            
        end
    end
end
%end

function parsave(fname, p,Y_t,t,simTime,y0)
save(fname, 'p', 'Y_t','t','simTime','y0');
end
