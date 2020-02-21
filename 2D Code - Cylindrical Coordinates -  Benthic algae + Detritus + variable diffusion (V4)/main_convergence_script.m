
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
parfor p = 1:3
    
    
    tic
    %% Model Parameters
    % all parameters are organized in a struct p for implementation convenience
    pp = struct;
    %%% Diffusion coefficients are defined below %%%
    pp.I0 = 300;    % Light intensity at the surface [micro-mol photons m^-2 s^-1]
    pp.kA = 0.0003; % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
    pp.kD = 0.0003; % Specific light-attenuation coefficient of detritus [m^2 mg C^-1]
    pp.kB = 0.00003; % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
    pp.kbg = 0.2;    % Background light-attenuation coefficient [m^-1]
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
    pp.r = 0.02;          % Specific mineralization rate of sedimented nutrients [day^-1]
    pp.vA = 0.1;          % Algal sinking speed [m day^-1]
    pp.vD = 0.25;         % detritus sinking speed [m day^-1]
    pp.Dbg = 0.02;        % remineralization of detritus in the water column [day^-1]
    pp.resus = 0.2;       % resuspension rate of the detritus from the sediment. [day^-1]
    pp.resus_Max = 0.02;  % Maximum resuspension rate of the detritus in the sediment [day^-1]
    pp.resus_H = 8;       % half saturation constant of the resuspension rate functional response [m^2 day^-1]
    pp.death_rate = 1;    % coefficient governing the proportion of sinking algae at the bottom that dies.
    % 0 = no death. 1 = all algae that would have sunk through the sediment dies.
    pp.benth_recycling = 0.5; % range: [0,1]. Governs the portion of respired nutrients that are released as dissolved nutrients.
    % the rest is bound in particulate matter in the sediment.
    
    
    
    %% Lake topology and Mesh
    % Quantities relating to system size
    res_vec = [4,16,32];
    pp.Xn   = res_vec(p)+1;  % Number of grid-points (width)
    pp.Zn   = res_vec(p)+1;  % Number of grid-points (depth)
    pp.Lmin = 0.1; % Minimum lake depth (depth at land-water interface) [m]
    pp.Lmax = 20;  % Maximum lake depth [m]
    pp.W    = 20;  % Lake radius [m]
    pp.alpha = 1;  % Exponent governing the slope of the lake bottom
    
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
    for i=1:pp.Zn
        for j = 1:pp.Xn
            pp.Z(i,j) = (pp.Lmax/(pp.Zn-1))*(i-1)*(1 + (pp.Lmin/pp.Lmax -1)* (pp.X(i,j)/pp.W).^(pp.alpha));
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
    
    pp.volumes_cyl = vol_areas_cyl_3d_fn(pp);
    pp.Area_bottom_cyl = Area_bottom_cyl_fn(pp);
    
    
    
    %% Diffusion Coefficients
    dx = zeros(pp.Zn-1,pp.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
    dz = zeros(pp.Zn-1,pp.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
    
    diff_max_x = 100; % horizontal diffusion coefficient at the surface
    diff_max_z = 100; % vertical diffusion coefficient at the surface
    diff_min_x = diff_max_x % horizontal diffusion coefficient at the bottom
    diff_min_z = diff_max_z % vertical diffusion coefficient at the bottom
    
    
    % Linearly decreasing diffusion coefficients with the lake depth.
    for j = 1:pp.Xn-1
        for i = 1:pp.Zn-1
            dx(i,j) =  diff_min_x + (diff_max_x-diff_min_x)*(1- pp.Z_vol(i,j)/pp.Lmax);
            dz(i,j) =  diff_min_z + (diff_max_z-diff_min_z)*(1- pp.Z_vol(i,j)/pp.Lmax);
        end
    end
    
    pp.dx = dx;
    pp.dz = dz;    
    
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
    pp.ntot_0 = ntot_0
    %% state variables
    A = A0'; % transpose of matrices in order to use colon notation to reshape to vector form.
    Rd = Rd0';
    D = D0';
    Rs = Rs0';
    B = B0';
    y0 = [A(:); Rd(:); D(:); Rs(:); B(:)];
    
    %% Simulation of model
    tend = 10; % end simulation time
    %t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
    %t_span = [1:tend];
    %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events',@(t,y) eventfun_V4(t,y,p) , 'NonNegative',(1:length(y0)));
    ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'NonNegative',find(y0));
    %ode_opts = odeset( 'abstol' , 1e-9 , 'reltol' , 1e-9 , 'Events', @(t,y) eventfun_V4(t,y,pp), 'NonNegative',(1:length(y0)));
    [t,Y_t] = ode15s( @(t,Y) rhs_function_V4(t,Y,pp), [0, tend] , y0 , ode_opts);
    %[t,Y_t] = ode15s( @(t,Y) dAdt_efficient_correct_V4_detritus(t,Y,p) , t_span , y0 , ode_opts);
    
    %% recording of simulation time & saving workspace
    simTime = toc;
    file_name = "2D_convergence_res_V4_dx_" + "Xn_" + num2str(pp.Xn) + "_Zn_" + num2str(pp.Zn) + "_end_time_" + num2str(tend);
    parsave(file_name,pp,Y_t,t,simTime,y0);
    
end

function parsave(fname, p,Y_t,t,simTime,y0)
save(fname, 'p', 'Y_t','t','simTime','y0');
end
