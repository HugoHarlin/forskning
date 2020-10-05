
%% a clean slate
%close all
%clear
%clc
tic

%% Implementation details

% Version 6 - detritus + variable diffusion coefficients. square grid, no
% transform! pelagic algae and benthic algae are one and the same species.
% When pelagic algae sinks to the bottom it doesn't die like in previous
% version but is instead "turned" into benthic algae, think of it as algae
% simply accumulating on the bottom and continuing to live on.

% This code implements a system of equations based on the work by Jäger et
% al. (2010) with an additional benthic algal layer at the bottom,
% in two spatial dimensions, using a finite volume method and Gauss divergence
% theorem to solve the system of equations.
% ode15s is used to iterate in time.
%
% Plankton nutrient density stochiometry q is constant, and the lake topography
% is described by an upsidedown cone where the bot37tom edge is sloped so
% that the lake is deepest in the middle and shallow at the shore. The slope of the
% bottom is governed by an exponent alpha, and half the lake is modelled so that the
% deepest edge corresponds the middle of the lake.
%
% The matrices in this script are indexed so that origo is centered at the center of the lake surface,
% with the z-axis facing downwards and the x-axis facing the lake shore.
%
% Written by Hugo Harlin 2020


%% Loop over different turbidities
parfor p =  1:3
    
    %% all parameters are organized in a struct pp for implementation convenience
    pp = struct;
    
    pp.perc_increase = p;
    resus_index = 3;
    kbg_index = 2;
    
    %for manual_resus_indx = 1 % :3
    
    manual_resus_vec =  [0.001, 0.01, 0.1];
    
    depth_var =  1; % thermocline depth index
    %p = 1;
    
    % for kbg_index = 1:28 %:2 % :3
    %for depth_var =  1:3  % looping over thermocline depth
    %for resus_index = 5 %2 :4
    %depth_var = 1;
    % for x_res_index = [0,1]
    %for therm_depth = [2,5,10,15]
    % kbg_index =  p;
    %kbg_index = 1;
    tic
    
    
    
    pp.manual_therm_resus_val = manual_resus_vec(manual_resus_indx);
    
    %% Lake topology and Mesh
    % Quantities relating to system size
    pp.Xn   = 21;  % Number of grid-points (width)
    pp.Zn   = 21;  % Number of grid-points (depth)
    
    %   res_vec = (10:4:60);
    % pp.Xn   = res_vec(p);  % Number of grid-points (width)
    %  pp.Zn   = res_vec(p);  % Number of grid-points (depth)
    
    pp.Lmin = 0; % Minimum lake depth (depth at land-water interface) [m]
    depth_vec = [1,2,5, 10, 20 40];
    %depth_vec = [5, 100];
    pp.Lmax = depth_vec(5); %20;  % Maximum lake depth [m]
    width_vec = [10,20,40,50,100,200,300,400];
    pp.W    =   width_vec(3);  % Lake radius [m]
    alpha_vec = [1, 1.5]; % slopes being looped over.
    alpha_index = 1;
    pp.alpha = alpha_vec(alpha_index);  % Exponent governing the slope of the lake bottom
    
    pp.stratified = 0; % if true, the lake diffusion coefficients are set so that the diffusion coefficient reflect a thermally stratified
    
    if(pp.stratified)
        pp.increased_x_res = 0;
        pp.increased_z_res = 0;
        % lake with an epilimnion at 5 meters, (3 meters wide transition)
        %depth_vec = [2, 5, 10, 15];
        depth_vec = [1, 4, 8, 10, 15];
        
        pp.thermocline_depth = depth_vec(depth_var); % vertical depth at which the water transitions from a well mixed state, the upper boundary of the thermocline [m].
        pp.thermocline_thickness = 3; % the thickness of the thermocline layer [m]
        pp.diff_above_thermocline = 100; % [day^-1]
        pp.diff_in_thermocline = 1; % [day^-1]
        pp.diff_below_thermocline = 10;   % [day^-1]
    end
    
    % Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
    % at the center of the lake (slope = alpha* (Lmin - Lmax)/W).
    % (0,0) is placed at the center of the lake at the surface, y-dim is facing
    % downward and x-dim is facing towards the lake edge
    
    pp.X = zeros(pp.Zn, pp.Xn); % Mesh-spacing in x-dimension
    pp.Z = zeros(pp.Zn, pp.Xn); % Mesh-spacing in y-dimension
    
    for ii=1:pp.Xn
        pp.X(ii,:) = (0:pp.W/(pp.Xn-1):pp.W)'; % even spacing horizontally
    end
    
    for ii=1:pp.Zn
        pp.Z(:,ii) = (0:pp.Lmax/(pp.Zn-1):pp.Lmax)';  % even spacing vertically
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
            pp.X =  [pp.X; pp.X(1:3* (pp.thermocline_thickness+1),:)];
            pp.Zn = pp.Zn + 3* (pp.thermocline_thickness+1); % updating the count of number of rows in our grid p.Zn
        end
        
        %             if(pp.increased_x_res) % horizontal grid refinement (more columns) where the thermocline meets the bottom.
        %                 x_strat_top = pp.W*((pp.thermocline_depth - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the top of the thermocline intersects with the bottom.
        %                 x_strat_bottom = pp.W*((pp.thermocline_depth + pp.thermocline_thickness - pp.Lmax)/(pp.Lmin - pp.Lmax))^(1/pp.alpha); % this is the depth at which the bottom of the thermocline intersects with the bottom.
        %                 x_strat_top = ceil(x_strat_top); % by rounding up to the nearest whole number, ensuring that the variable is compatible as an index and that we dont miss any grid refinement
        %                 x_strat_bottom = floor(x_strat_bottom); % where the thermocline is.
        %
        %
        %                 for ii=1:1:pp.Zn
        %                     pp.X(ii,:) = [0:pp.W/(pp.Xn-1):pp.W]; % even spacing of the grid in x-dimension
        %                 end
        %
        %                 % the grid is compressed in y-dimension, with depth Lmin at the
        %                 % shore and Lmax at the center of the lake.
        %                 tt=1:pp.Zn;
        %
        %                 temp_index = x_strat_bottom +1 ; % index used for concatenating existing Z-matrix with new rows
        %
        %                 for jjj = x_strat_bottom:x_strat_top-1
        %                     %for j = [1]
        %                     num_new_columns = 1;
        %                     new_columns = zeros(pp.Zn,num_new_columns); % new columnds for pp.Z matrix
        %                     temp_X = zeros(pp.Zn,num_new_columns); % new colums for pp.X matrix
        %
        %                     for k = 1:num_new_columns
        %                         new_columns(:,k) = (tt-1)'.*(pp.Lmax/(pp.Zn-1)).*(1 + (pp.Lmin/pp.Lmax -1).*((pp.X(tt, temp_index)+k/( num_new_columns +1) )/pp.W).^(pp.alpha));
        %                         temp_X(:,k) = pp.X(tt,temp_index) + k/(num_new_columns+1);
        %                     end
        %
        %                     pp.X  = [pp.X(:,1:temp_index)  temp_X  pp.X(:,temp_index+1:end)];
        %
        %                     pp.Z = [pp.Z(:,1:temp_index) new_columns pp.Z(:,temp_index+1:end)]; % adding in the new columns in the Z-matrix
        %                     temp_index = temp_index +num_new_columns+1; % adjusting the index where the new rows are to be added.s
        %                 end
        %
        %                 idx = jjj-x_strat_bottom+1;
        %                 pp.Xn = pp.Xn + num_new_columns* (idx);
        %             end
        
    end
    
    % coordinates of the center of each grid element
    X_vol = zeros(pp.Zn-1, pp.Xn-1);
    Z_vol = zeros(pp.Zn-1, pp.Xn-1);
    
    for jjj=1:pp.Xn-1
        for ii =1:pp.Zn-1
            x_idx = jjj+1; % Required (static): When using the nested for-loop variable for indexing a sliced array, you must use the variable in plain form, not as part of an expression.
            z_idx = ii+1;
            X_vol(ii,jjj) = (pp.X(ii,jjj) + pp.X(z_idx,jjj) + pp.X(ii,x_idx) + pp.X(z_idx,x_idx) ) /4;
            Z_vol(ii,jjj) = (pp.Z(ii,jjj) + pp.Z(z_idx,jjj) + pp.Z(ii,x_idx) + pp.Z(z_idx,x_idx) ) /4;
        end
    end
    
    
    % removing nodes outside the grid.
    jjjj = pp.Xn-1;
    for ii=2:pp.Zn-1
        Z_vol(ii,jjjj:end) = NaN;
        X_vol(ii,jjjj:end) = NaN;
        jjjj = jjjj -1;
    end
    
    % Saving values in the struct
    pp.X_vol = X_vol;
    pp.Z_vol = Z_vol;
    
    
    % Calculating the volume of each grid element, assuming rotational symmetry
    % around the z-axis through the lake enter.
    
    volumes_cyl = zeros(pp.Zn-1,pp.Xn-1);
    
    for ii =1:pp.Zn-1
        for jjj =1:pp.Xn-1
            
            x_idx = jjj+1; % loop index is not allowed as part of an expression when indexing array in parfor (sigh)
            z_idx = ii+1;
            
            x1 = pp.X(ii,jjj);
            x2 = pp.X(z_idx,jjj);
            x3 = pp.X(z_idx,x_idx);
            x4 = pp.X(ii,x_idx);
            
            Z1 = pp.Z(ii,jjj);
            Z2 = pp.Z(z_idx,jjj);
            Z3 = pp.Z(z_idx,x_idx);
            Z4 = pp.Z(ii,x_idx);
            
            volumes_cyl(ii,jjj) = 2*pi*( ((x4^3)/3 - (x1^3)/3)*(Z3-Z2-Z4+Z1)/(x3-x2) + ((x4^2)/2 - (x2^2)/2)*(Z2-Z1 - x2*(Z3-Z2-Z4+Z1)/(x3-x2)));
        end
    end
    
    jj = pp.Xn-1;
    for ii=2:pp.Zn-1
        volumes_cyl(ii,jj:end) = NaN;  % removing nodes outside the grid.
        jj = jj -1;
    end
    
    jj = pp.Xn-1;
    for ii =1:pp.Zn-1
        
        z_idx = ii+1;
        x_idx = jj+1;
        
        x1 = pp.X(ii,jj);   % top left
        x2 = pp.X(z_idx,jj); % bottom
        x4 = pp.X(ii,x_idx); % top right
        
        Z1 = pp.Z(ii,jj);   % top left
        Z2 = pp.Z(z_idx,jj); % bottom
        Z4 = pp.Z(ii,x_idx); % top right
        
        slope = abs(Z4-Z2)/(x4-x2); % slope of the bottom segment
        
        volumes_cyl(ii,jj) =   pi/(slope^2)*(1/3*(Z2-Z4)^3 - (pp.Lmax-Z1)*(Z2-Z4)^2 + (pp.Lmax-Z1)^2*(Z2-Z4));% Diagonal triangular elements rotated around the z-axis
        r_circ = (Z2-Z4)*pi*x1^2;
        volumes_cyl(ii,jj) = volumes_cyl(ii,jj) - r_circ;
        
        jj = jj-1;
    end
    pp.volumes_cyl = volumes_cyl(~isnan(volumes_cyl));
    pp.Area_bottom_cyl = Area_bottom_cyl_fn(pp); % area of each grid element on the bottom
    
    %% Model Parameters
    
    pp.model_version = 6; % choose between 5 and 6
    % currently (sep 2020) model 5 and 6 are implemented.
    
    % model 5: benthic algae and pelagic algae different
    % species. Algae die instantly when it sinks into the bottom,
    % and is turned into sedimented nutrients. Benthic algae does
    % not resuspend, only sediment as detritus.
    
    % model 6: Benthic algae and pelagic algae are one and the same
    % species, when pelagic algae sinks into the bottom it is
    % turned into benthic algae. Benthic algae resuspends, and is
    % turned into pelagic algae when it happens.
    
    
    %%% Diffusion coefficients are defined below %%%
    pp.I0 = 300;    % Light intensity at the surface [micro-mol photons m^-2 s^-1]
    pp.kA = 0.0003; % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
    pp.kD = 0.00003; % Specific light-attenuation coefficient of detritus [m^2 mg C^-1]
    pp.kB = 0.0003; % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
    %kgb_vec = [ 0.4 0.6 2.0];
    kgb_vec = [ 0.2 0.6 0.8];
    %kgb_vec = [0 0.1 0.2 0.35 0.6 1.0 2.0];
    pp.kbg = kgb_vec(kbg_index);  % Background light-attenuation coefficient [m^-1]
    pp.lbg_A =  0.08; % Specific algal maintenance respiration losses [day^-1]
    pp.Ad = 0.02;    % algal death rate [day^-1]
    pp.M = 1.5;      % Half-saturation constant of algal nutrient uptake [mg P m^-3]
    pp.M_benth = 1.5; % Half-saturation constant of benthic algae nutrient uptake [mgPm^-3]
    pp.Gmax =  1.08;   % Maximum specific phytoplankton production rate [day^-1]
    pp.Gmax_benth = 1.08; % Maximum specific benthic algae production rate [day^-1]
    pp.H = 120.0;         % Half-saturation constant of light-dependent algal production [micro-mol photons m^-2 s^-1]
    pp.H_benth = 120.0;   % Half-saturation constant of light-dependent benthic algal production [micro-mol photons m^-2 s^-1]
    pp.q = 0.0244;        % Algal nutrient quota, Redfield ratio [mgP/mgC]
    pp.lbg_benth = 0.1;   % Specific benthic algae maintenance respiration losses [day^-1]
    remin_vec = [0, 0.01, 0.02, 0.05, 0.1, 0.2]; % 0.02 is used as default value
    remin_index = 3; %3;
    pp.r = remin_vec(remin_index);  % Specific mineralization rate of sedimented nutrients [day^-1]
    pp.vA = 0.1;          % Algal sinking speed [m day^-1]
    pp.vD = 0.25;     % detritus sinking speed [m day^-1]
    %defined below.  pp.resus_benth = 0.2; % resuspension rate of the benthic algae, resuspended benthic algae becomes pelagic algae(its treated as one and the same species)
    pp.Dbg =  0.02;       % remineralization of detritus in the water column [day^-1]
    pp.resus_Max = 0.02;  % Maximum resuspension rate of the detritus in the sediment [day^-1]
    pp.resus_H = 8;       % half saturation constant of the resuspension rate functional response [m^2 day^-1]
    % 0 = no death. 1 = all algae that would have sunk through the sediment dies.
    pp.benth_recycling = 0.8; % 0.8 used as standard. range: [0,1]. Governs the portion of respired nutrients that are released as dissolved nutrients.
    % the rest is bound in particulate matter in the sediment.
    
    
    %%% Seasonality
    pp.seasonality = 0; % if true the simulation is run seasonally untill a steady state is reached or the simulation has run for 10000 years.
    pp.seasonality_light = 1; % seasonally varying light (sinusoidal)
    pp.maxLight = 400; % maximum light intesity for seasonally varying light [micro-mol photons m^-2 s^-1]
    pp.minLight = 200; % minimum light intesity for seasonally varying light [micro-mol photons m^-2 s^-1]
    pp.seasonality_thermoC = 0; % seasonally varying depth of the thermocline;
    pp.minTherm = 2 ; % shallowest depth of the thermocline (mid winter?) [m]
    pp.maxTherm = 10; % greatest depth of the thermocline (mid summer?) [m]
    pp.seasonality_mixing = 0; %set to true if periodic mixing is desired (discrete mixing event)
    pp.mixing_frequency = 365; % time between each mixing event [days]
    pp.max_sim_time = 365*1e+4; % upper bound of the total simulation time if steady state is not reached before then (10000 years).
    
    %% looping over the diffusion coefficients
    % for diff_index = 1:3 %,1000]
    %for diff_max_z = [1,10,100] %,1000]
    
    %% Diffusion Coefficients
    dx = zeros(pp.Zn-1,pp.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
    dz = zeros(pp.Zn-1,pp.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
    
    if(~pp.stratified)
        % diff_vec =[0.001, 0.1, 1, 5, 10,100, 1000];
        diff_vec =[0.1, 1, 5, 10,100, 1000];
        %diff_max_x =diff_vec(p);
        diff_index = 4;
        diff_max_x = diff_vec(4); % horizontal diffusion coefficient at the surface
        %diff_max_z = 100; % vertical diffusion coefficient at the surface
        diff_min_x =  diff_vec(4); % horizontal diffusion coefficient at the bottom
        diff_min_z =  diff_vec(diff_index); % diff_max_x; % vertical diffusion coefficient at the bottom
        diff_max_z = diff_min_z; % diff_max_x;
        %diff_min_z = diff_max_z; % vertical diffusion coefficient at the bottom
        
        dx = diff_max_x.*ones(pp.Zn-1,pp.Xn-1);
        dz = diff_max_z.*ones(pp.Zn-1,pp.Xn-1);
        
        
        % Linearly decreasing diffusion coefficients with the lake depth.
        %                 for j = 1:pp.Xn-1
        %                     for ty = 1:pp.Zn-1
        %                         dx(ty,j) =  diff_min_x + (diff_max_x-diff_min_x)*(1- pp.Z_vol(ty,j)/pp.Lmax);
        %                         dz(ty,j) =  diff_min_z + (diff_max_z-diff_min_z)*(1- pp.Z_vol(ty,j)/pp.Lmax);
        %                     end
        %                 end
        
        pp.dx = dx;
        pp.dz = dz;
    end
    
    
    
    if(pp.stratified)
        pp.dx = zeros(pp.Zn-1);
        pp.dz = zeros(pp.Xn-1);
        %%%%% manual setting of the diffusion coefficients %%%%
        for x = 1:pp.Xn-1
            for z = 1:pp.Zn-1
                
                if(pp.Z_vol(z,x) < pp.thermocline_depth)
                    
                    % pp.dx(z,x) = pp.diff_above_thermocline;
                    % pp.dz(z,x) =  pp.diff_above_thermocline; % *(pp.Zn-z)/pp.Zn;
                    
                    % linearly decreasing with depth to 80% of the max
                    % value
                    pp.dx(z,x) =  pp.diff_above_thermocline *(0.8 + 0.2*( pp.thermocline_depth- abs(pp.Z_vol(z,x)))/pp.thermocline_depth); %
                    pp.dz(z,x) =  pp.diff_above_thermocline *(0.8 + 0.2*( pp.thermocline_depth- abs(pp.Z_vol(z,x)))/pp.thermocline_depth); %
                    
                elseif(pp.thermocline_depth <= pp.Z_vol(z,x) &&  pp.Z(z,x) < pp.thermocline_depth + pp.thermocline_thickness)
                    pp.dx(z,x) = pp.diff_in_thermocline;
                    pp.dz(z,x) = pp.diff_in_thermocline;
                else
                    
                    pp.dx(z,x) =  pp.diff_below_thermocline; % *(pp.Zn-0.5*z)/(pp.Zn);
                    pp.dz(z,x) =  pp.diff_below_thermocline; % *(pp.Zn-0.5*z)/(pp.Zn);
                    
                    
                    pp.dx(z,x) =  pp.diff_below_thermocline *(0.8 + 0.2*( pp.Lmax- abs(pp.Z_vol(z,x)))/pp.Lmax) ; %
                    pp.dz(z,x) =  pp.diff_below_thermocline *(0.8 + 0.2*( pp.Lmax- abs(pp.Z_vol(z,x)))/pp.Lmax) ; %
                    
                end
            end
        end
    end
    
    
    
    
    % removing values outside the grid
    
    jj = pp.Xn-1;
    for ii=2:pp.Zn-1
        pp.dx(ii,jj:end) = NaN;
        pp.dz(ii,jj:end) = NaN;
        jj = jj -1;
    end
    
    
    test = 1;
    
    %% resuspension
    
    pp.constant_resuspension = 1;
    
    pp.manual_in_therm_resus = false;
    % pp.response_type =  p;
    pp.response_type = 1;
    pp.stepFun_resus = false; % true;
    
    
    if(pp.constant_resuspension)
        resus_vec = [0 0.01 0.05, 0.1, 0.2, 0.3 0.6]; % 0.05 is used as default value
        % resus_vec = [0.01 0.02 0.03 .04] % 0.05 is used as default value
        % resus_vec = [0.001 0.003 0.006] % 0.05 is used as default value
        % resus_vec = [0.001:0.005:0.1 0.15 : 0.05 : 0.5];
        %resus_index = 3;
        pp.resus = ones(1,pp.Xn-1).*resus_vec(resus_index);  % resuspension rate of the detritus from the sediment. [day^-1]
        pp.resus_benth = pp.resus;  % same resuspension rate for detritus and benthic algae
    else
        
        
        %%%% Resuspension rate that varies with depth %%%%
        pp.resus = ones(1,pp.Xn-1);
        
        % Type 1 functional response (linear)
        if(pp.response_type == 1)
            max_resus = 0.2;
            min_resus = 0;
            b_idx = pp.Zn-1; % bottom index
            for i=1:pp.Xn-1
                %   diff_normal = (pp.dx(i,b_idx) + pp.dz(i,b_idx)* pp.W/pp.Lmax)/sqrt(1 + pp.W^2/pp.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                diff_normal = pp.dz(b_idx,i);
                pp.resus(i) = min_resus +  ((max_resus - min_resus)/100)*diff_normal;
                b_idx = pp.Zn-1 -i;
            end
        end
        
        % Type 2 functional response
        if(pp.response_type == 2)
            % the parameters are tweaked manually to get a
            % functional response that had a "good" s-shape to
            % it, assuming a maimum resuspension rate of 20%/day.
            
            c =1; % 4.2;
            b = 0.3; %2.0e-05;
            
            b_idx = pp.Zn-1; % bottom index
            for i=1:pp.Xn-1
                %   diff_normal = (pp.dx(i,b_idx) + pp.dz(i,b_idx)* pp.W/pp.Lmax)/sqrt(1 + pp.W^2/pp.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                diff_normal = pp.dz(b_idx,i);
                pp.resus(i) = (0.2*b.*diff_normal.^c)./(1 + 1*b.*diff_normal.^c);
                b_idx = pp.Zn-1 -i;
            end
            test = 1;
        end
        
        
        % Type 3 functional response
        if(pp.response_type == 3)
            % the parameters are tweaked manually to get a
            % functional response that had a "good" s-shape to
            % it, assuming a maimum resuspension rate of 20%/day.
            
            c =3; % 4.2;
            b = 0.02; %2.0e-05;
            
            b_idx = pp.Zn-1; % bottom index
            for i=1:pp.Xn-1
                %   diff_normal = (pp.dx(i,b_idx) + pp.dz(i,b_idx)* pp.W/pp.Lmax)/sqrt(1 + pp.W^2/pp.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                diff_normal = pp.dz(b_idx,i);
                pp.resus(i) = (0.2*b.*diff_normal.^c)./(100 + 1*b.*diff_normal.^c);
                b_idx = pp.Zn-1 -i;
            end
        end
        
        
        if(pp.manual_in_therm_resus) % sets the resuspension coefficient in the thermocline manually
            b_idx = pp.Zn-1; % bottom index
            for i=1:pp.Xn-1
                %   diff_normal = (pp.dx(i,b_idx) + pp.dz(i,b_idx)* pp.W/pp.Lmax)/sqrt(1 + pp.W^2/pp.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                if(pp.dz(b_idx,i) == pp.diff_in_thermocline)
                    pp.resus(i) = pp.manual_therm_resus_val;
                end
                b_idx = pp.Zn-1 -i;
            end
            
        end
        
        
        % stepfunction resuspension
        if(pp.stepFun_resus)
            
            pp.upper_resusp = 0.05; % resuspension rate above upper threshold [1/day]
            pp.lower_resusp = 0.0001; % resuspension rate below lower threshold [1/day]
            
            pp.upper_thresh = 10; % upper threshold depth [m]
            pp.lower_thresh = 12; % lower threshold depth [m]
            
            
            for i=1:pp.Xn-1
                if(pp.Z_vol(pp.Xn-i) < pp.upper_thresh) % above upper threshold.
                    pp.resus(i) = pp.upper_resusp;
                    
                    
                elseif(pp.Z_vol(pp.Xn-i) > pp.lower_thresh) % below lower threshold.
                    pp.resus(i) = pp.lower_resusp;
                    
                else % we are between the thresholds
                    
                    % diff_normal = (pp.dx(i,b_idx) + pp.dz(i,b_idx)* pp.W/pp.Lmax)/sqrt(1 + pp.W^2/pp.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                    % diff_normal = pp.dz(b_idx,i);
                    pp.resus(i) =   pp.lower_resusp +    (pp.upper_resusp - pp.lower_resusp)*(1- ( pp.Z_vol(pp.Xn-i)- pp.upper_thresh )/(pp.lower_thresh - pp.upper_thresh));
                    
                    
                end
            end
        end
        
        pp.resus_benth = pp.resus;
    end
    
    %% Inital Conditions
    % A0  = 1.0*ones(pp.Zn-1, pp.Xn-1);    % Initial Algal carbon density [mg C m^-3]
    % Rd0 = 10.000*ones(pp.Zn-1, pp.Xn-1); % initial concentration of dissolved nutrients [mg P m^-3]
    % D0  = 1.000*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of detritus [mg P m^-3]
    % Rs0 = 1.0*ones(1, pp.Xn-1);          % initial concentration of sediment nutrient density [mg P m^-2]
    % B0  = 1.00*ones(1, pp.Xn-1);         % initial concentration of benthic algal density [mg C m^-2]
    
    % constant ammount of nutrients per surface area
    nutr_area      = (9.4178e+04)./(pi*20^2); % [mgP/m^2] corresponds to a lake with radius 20m, depth 20m, A0 = 1 mgC/(m^3), Rd0 = 10 mgP/(m^3), D0 =1 mgP/(m^3), Rs0 = 1 mgP/(m^2), B0 = 1 mgP/(m^2)
     perc_vec = [0.2 1 5];
    nutr_tot_start =  perc_vec(pp.perc_increase)*(pi*pp.W^2)*nutr_area; % total amount of starting nutrients in lake [mgP/m^3]
    
    %             % percentages of total nutrients bound in each state variable
    perc_A0  = 204.4130/9.4178e+04; % percentage of total starting nutrients in pelagic algae
    perc_Rd0 = 8.3776e+04/9.4178e+04; % percentage of total starting nutrients in pelagic dissolved nutrientss
    perc_D0  = 8.3776e+03/9.4178e+04; % percentage of total starting nutrients in pelagic detritus
    perc_B0  = 43.3625/9.4178e+04; % percentage of total starting nutrients in benthic algae
    perc_Rs0 = 1.7772e+03/9.4178e+04; % percentage of total starting nutrients in sedimented nutrients
    
    
    %             perc_A0  = 0.66; % percentage of total starting nutrients in pelagic algae
    %             perc_Rd0 = 0.31; % percentage of total starting nutrients in pelagic dissolved nutrientss
    %             perc_D0  = 0.01; % percentage of total starting nutrients in pelagic detritus
    %             perc_B0  = 0.01; % percentage of total starting nutrients in benthic algae
    %             perc_Rs0 = 0.01; % percentage of total starting nutrients in sedimented nutrients
    %
    %
    
   
%     perc_A0  =   perc_A0; % percentage of total starting nutrients in pelagic algae
%     perc_Rd0 =   perc_Rd0; % percentage of total starting nutrients in pelagic dissolved nutrientss
%     perc_D0  =   perc_D0; % percentage of total starting nutrients in pelagic detritus
%     perc_B0  =   perc_B0; % percentage of total starting nutrients in benthic algae
%     perc_Rs0 =   perc_Rs0; % percentage of total starting nutrients in sedimented nutrients
%     
%     
%     
    
    % concentratrions of state variables
    A0  = perc_A0/pp.q*nutr_tot_start/sum(pp.volumes_cyl);     % [mgC/m^3]
    Rd0 = perc_Rd0*nutr_tot_start/sum(pp.volumes_cyl);         % [mgP/m^3]
    D0  = perc_D0*nutr_tot_start/sum(pp.volumes_cyl);          % [mgP/m^3]
    B0  = perc_B0/pp.q*nutr_tot_start/sum(pp.Area_bottom_cyl); % [mgP/m^2]
    Rs0 = perc_Rs0*nutr_tot_start/sum(pp.Area_bottom_cyl);     % [mgP/m^2]
    
    %
    A0  = A0*ones(pp.Zn-1, pp.Xn-1);    % Initial Algal carbon density [mg C m^-3]
    Rd0 = Rd0*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of dissolved nutrients [mg P m^-3]
    D0  = D0*ones(pp.Zn-1, pp.Xn-1);  % initial concentration of detritus [mg P m^-3]
    B0  = B0*ones(1, pp.Xn-1);        % initial concentration of benthic algal density [mg C m^-2]
    Rs0 = Rs0*ones(1, pp.Xn-1);        % initial concentration of sediment nutrient density [mg P m^-2]
    
    
    
    % removing values outside grid
    jj= pp.Xn-1;
    for ii=2:pp.Zn-1
        A0(ii,jj:end) = NaN; % removing nodes outside the grid.
        Rd0(ii,jj:end) = NaN; % removing nodes outside the grid.
        D0(ii,jj:end) = NaN; % removing nodes outside the grid.
        jj = jj -1;
    end
    
    A0 = A0(~isnan(A0));
    Rd0 = Rd0(~isnan(Rd0));
    D0 = D0(~isnan(D0));
    
    %% calculation of total nutrient content at t=0, (to be conserved)
    n_algae_0 = pp.volumes_cyl.*A0*pp.q;
    n_dissolved_0 = pp.volumes_cyl.*Rd0;
    n_detritus_0 = pp.volumes_cyl.*D0;
    n_sediment_0 = Rs0.*pp.Area_bottom_cyl;
    n_benthic_0 = B0.*pp.q.*pp.Area_bottom_cyl;
    
    pp.ntot_algae_0  = sum(sum(n_algae_0));
    pp.ntot_dissolved_0 = sum(sum(n_dissolved_0));
    pp.ntot_detritus_0 = sum(sum(n_detritus_0));
    pp.ntot_sed_0 = sum(n_sediment_0);
    pp.ntot_benthic_0 = sum(n_benthic_0);
    
    ntot_0 =  pp.ntot_algae_0 + pp.ntot_dissolved_0 + pp.ntot_detritus_0 + pp.ntot_sed_0 + pp.ntot_benthic_0;
    pp.ntot_0 =  ntot_0;
    
    %% state variables
    A = A0; % transpose of matrices in order to use colon notation to reshape to vector form.
    Rd = Rd0;
    D = D0;
    Rs = Rs0;
    B = B0;
    y0 = [A(~isnan(A)); Rd(~isnan(Rd));D(~isnan(D)); Rs(:); B(:)];
    
    %% Simulation of model
    
    
    ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events', @(t,y) eventfun_V6(t,y,pp), 'NonNegative',(1:length(y0)));
    
    
    if(~pp.seasonality) % no periodic mixing
        %tend = 2000; % end simulation time
        %t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
        %t_span = [1:tend];
        %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events',@(t,y) eventfun_V6(t,y,p) , 'NonNegative',(1:length(y0)));
        %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'NonNegative',find(y0));
        [t,Y_t] = ode15s( @(t,Y) rhs_function_V6(t,Y,pp), [0 inf] , y0 , ode_opts);
        %[t,Y_t] = ode15s( @(t,Y) rhs_function_V6(t,Y,pp), [0 10000] , y0 , ode_opts);
        %[t,Y_t] = ode15s( @(t,Y) dAdt_efficient_correct_V6_detritus(t,Y,p) , t_span , y0 , ode_opts);
        
    else % Seasonality active
        t_var = 0; % variable keeping track of the total lapsed simulation time
        res_total = []; % state variables for all time states of the total simulation, analog of Y_t for a single ode15s run
        time_total = []; % this vector will contain the simulation time for the total simulation, analog to the t vector if ode15s is only run once.
        
        while(t_var <  pp.max_sim_time)
            [t,Y_t] = ode15s( @(t,Y) rhs_function_V6(t,Y,pp), [ t_var:1:(t_var + pp.mixing_frequency)] , y0 , ode_opts);
            
            % saving timesteps and state variables
            time_total = [time_total; t];
            res_total = [res_total; Y_t];
            
            % updating the timekeeping variable
            t_var = t_var +pp.mixing_frequency;
            
            [xx,isterm,dir] = eventfun_V6(t(end),Y_t(end,:)',pp);
            
            if(xx >0) % steady state is not reached and we have arrived at the mixing event
                %%%%%%%%%%% mixing event %%%%%%%%%%%%
                % Extracting results
                
                A = Y_t(end,1:sum(1:pp.Xn-1)); % phytoplankton
                Rd = Y_t(end,sum(1:pp.Xn-1)+1:2*sum(1:pp.Xn-1)); % dissolved nutrients
                D = Y_t(end,2*sum(1:pp.Xn-1) +1 : 3*sum(1:pp.Xn-1)); % ditritus
                Rs = Y_t(end,3*sum(1:pp.Xn-1)+1 : 3*sum(1:pp.Xn-1)+ pp.Xn-1); % sedimented nutrients
                B = Y_t(end, 3*sum(1:pp.Xn-1)+ pp.Xn :  3*sum(1:pp.Xn-1)+ 2*(pp.Xn-1)); % benthic algae
                
                
                % if seasonal mixing is set to true the code below
                % mixes the lake.
                if(pp.seasonality_mixing)
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
                end
                
                y0 = [A(:); Rd(:); D(:); Rs(:); B(:)]; % state variable input for the next function call of ode15s, with the bulk of the lake being well mixed.
                
            else
                t_var = pp.max_sim_time +1; % we exit the loop if steady state is reached.
            end
            
            
            
            y1 = (res_total(end, :));
            y2 = (res_total(end-365, :));
            
            if(norm(y2-y1) < 1e-6)
                t_var = pp.max_sim_time +1;  % we exit the loop because the results are cyclic
            end
        end
        Y_t = res_total;
        t = time_total;
        
        
        
    end
    
    %% extracting results and calculating ntot end
    A  = Y_t(end,1:sum(1:pp.Xn-1)); % phytoplankton
    Rd = Y_t(end,sum(1:pp.Xn-1)+1:2*sum(1:pp.Xn-1)); % dissolved nutrients
    D  = Y_t(end,2*sum(1:pp.Xn-1) +1 : 3*sum(1:pp.Xn-1)); % ditritus
    Rs = Y_t(end,3*sum(1:pp.Xn-1)+1 : 3*sum(1:pp.Xn-1)+ pp.Xn-1); % sedimented nutrients
    B  = Y_t(end, 3*sum(1:pp.Xn-1)+ pp.Xn :  3*sum(1:pp.Xn-1)+ 2*(pp.Xn-1)); % benthic algae
    
    %Calculating nutrient content at t = Tend
    n_algae_end     = pp.volumes_cyl.*A'*pp.q;
    n_dissolved_end = pp.volumes_cyl.*Rd';
    n_detritus_end  = pp.volumes_cyl.*D';
    n_sediment_end  = Rs.*pp.Area_bottom_cyl;
    n_benthic_end   = B.*pp.q.*pp.Area_bottom_cyl;
    ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);
    
    %% recording of simulation time & saving workspace
    pp.simTime = toc;
    
    %% saving results
    if(1)
        if(pp.model_version == 6)
            % pelagic = benthic algae, pelagic sinking into
            % sediment layer turns into benthic algae. Benthic
            % algae resuspend into pelagic algae.
            file_name = "2D_benthic_results_V6_alpha_" + strrep(num2str(pp.alpha),'.','_') + "_kbg_" + strrep(num2str(pp.kbg),'.','_');
        end
        
        if(pp.model_version == 5)
            % Pelagic algae and benthic algae different species.
            % Algae die instantly when they sink into the bottom,
            % and turn into sedimented nutrients. Benthic algae are
            % not resuspended, only benthic algae.
            file_name = "2D_benthic_results_V5_alpha_" + strrep(num2str(pp.alpha),'.','_') + "_kbg_" + strrep(num2str(pp.kbg),'.','_');
        end
        
        
        file_name = file_name + "_res_" + strrep(num2str(pp.Xn),'.','_');
        %  file_name = file_name + "_max_depth_" + strrep(num2str(pp.Lmax),'.','_');
        
        if(pp.stratified)
            file_name = file_name +  "_Stratified_therm_depth_" + strrep(num2str(pp.thermocline_depth),'.','_');
            file_name = file_name +  "_thickness_" + strrep(num2str(pp.thermocline_thickness),'.','_');
            file_name = file_name + "_profile_" + num2str(pp.diff_above_thermocline)+ "_" + num2str(pp.diff_in_thermocline) + "_" + num2str(pp.diff_below_thermocline) + "";
            
            if(pp.increased_x_res)
                file_name = file_name + "_increased_x_res";
            end
        else
            file_name = file_name  +"_dx_" + strrep(num2str(diff_max_x),'.','_') + "_dz_" + strrep(num2str(diff_max_z),'.','_');
        end
        
        
        if(pp.seasonality)
            
            file_name = file_name +  "_seasonality_" ;
            if(pp.seasonality_thermoC)
                file_name = file_name +  "_therm_" ;
            end
            
            if(pp.seasonality_light)
                file_name = file_name +  "_light_" ;
            end
            
        end
        
        
        if(pp.constant_resuspension == 1)
            file_name = file_name + "_resus_rate_" + strrep(num2str(pp.resus(1)),'.','_') ;
        elseif(pp.stepFun_resus)
            file_name = file_name +  "_stepFun_resus_upper_th_" + strrep(num2str(pp.upper_thresh),'.','_') + "_lower_th_" +strrep(num2str(pp.lower_thresh),'.','_') + "_lower_resus_" + strrep(num2str(pp.lower_resusp),'.','_')+ "_lower_resus_ " +strrep(num2str(pp.upper_resusp),'.','_') ;
            
        else
            file_name = file_name +  "_variable_resus_resp_fn_" + num2str(pp.response_type);
            
            if(pp.manual_in_therm_resus)
                file_name = file_name +  "manual_in_therm_resus_" + strrep(num2str(pp.manual_therm_resus_val),'.','_');
            end
        end
        
        file_name = file_name + "_width_" +  strrep(num2str(pp.W),'.','_');
        file_name = file_name + "_depth_" +  strrep(num2str(pp.Lmax),'.','_');
        file_name = file_name + "_nutr_mult" +       strrep(num2str(perc_vec(pp.perc_increase)),'.','_');
        parsave(file_name,pp,Y_t,t,pp.simTime,y0);
        
    end
    
    % end
    
    % end
end
%end

function parsave(fname, p,Y_t,t,simTime,y0)
save(fname, 'p', 'Y_t','t','simTime','y0');
end

