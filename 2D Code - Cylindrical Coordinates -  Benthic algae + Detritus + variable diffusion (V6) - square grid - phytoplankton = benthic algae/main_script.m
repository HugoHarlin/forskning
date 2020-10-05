
%% Implementation details and a clean slate
close all
clear
clc
tic

% Version 4 - detritus + variable diffusion coefficients.

% This code implements a system of equations based on the work by J�ger et
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

%% Model Parameters
% all parameters are organized in a struct p for implementation convenience
p = struct;
%%% Diffusion coefficients are defined below %%%
p.I0 = 300;      % Light intensity at the surface [micro-mol photons m^-2 s^-1]
p.kA = 0.0003;   % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
p.kD = 0.0003;   % Specific light-attenuation coefficient of detritus [m^2 mg C^-1]
p.kB = 0.0003;   % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
p.kbg = 0.2;     % Background light-attenuation coefficient [m^-1]
p.lbg_A = 0.08;  % 0.08; % Specific algal maintenance respiration losses [day^-1]
p.Ad = 0.02;     % algal death rate [day^-1]
p.M = 1.5;       % Half-saturation constant of algal nutrient uptake [mg P m^-3]
p.M_benth = 1.5; % Half-saturation constant of benthic algae nutrient uptake [mgPm^-3]
p.Gmax = 1.08;   % Maximum specific phytoplankton production rate [day^-1]
p.Gmax_benth = 1.08; % Maximum specific benthic algae production rate [day^-1]
p.H = 120.0;         % Half-saturation constant of light-dependent algal production [micro-mol photons m^-2 s^-1]
p.H_benth = 120.0;   % Half-saturation constant of light-dependent benthic algal production [micro-mol photons m^-2 s^-1]
p.q = 0.0244;        % Algal nutrient quota, Redfield ratio [mgP/mgC]
p.q_benth =  0.0244; % Benthic algae nutrient quota, Redfield ratio [mgP/mgC]
p.lbg_benth = 0.1;   % Specific benthic algae maintenance respiration losses [day^-1]
p.r = 0.02;          % Specific mineralization rate of sedimented nutrients [day^-1]
p.vA = 0.1;          % Algal sinking speed [m day^-1]
p.vD = 0.25;         % detritus sinking speed [m day^-1]
p.Dbg = 0.02;        % remineralization of detritus in the water column [day^-1]
p.resus = 0.2;       % resuspension rate of the detritus from the sediment. [day^-1]
p.resus_Max = 0.02;  % Maximum resuspension rate of the detritus in the sediment [day^-1]
p.resus_H = 8;       % half saturation constant of the resuspension rate functional response [m^2 day^-1]
p.death_rate = 1;    % coefficient governing the proportion of sinking algae at the bottom that dies.
% 0 = no death. 1 = all algae that would have sunk through the sediment dies.
p.benth_recycling = 0.5; % range: [0,1]. Governs the portion of respired nutrients that are released as dissolved nutrients.
% the rest is bound in particulate matter in the sediment.

%% Lake topology and Mesh
% Quantities relating to system size
p.Xn   = 5;   % Number of grid-points (width)
p.Zn   = 5;   % Number of grid-points (depth)
p.Lmin = 0.1;  % Minimum lake depth (depth at land-water interface) [m]
p.Lmax = 20;   % Maximum lake depth [m]
p.W    = 20;   % Lake radius [m]
p.alpha = 0.5; % Exponent governing the slope of the lake bottom

% Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
% at the center of the lake (slope = alpha* (Lmin - Lmax)/W).
% (0,0) is placed at the center of the lake at the surface, y-dim is facing
% downward and x-dim is facing towards the lake edge

p.X = zeros(p.Zn, p.Xn); % Mesh-spacing in x-dimension
p.Z = zeros(p.Zn, p.Xn); % Mesh-spacing in y-dimension

for i=1:1:p.Zn
    p.X(i,:) = [0:p.W/(p.Xn-1):p.W]; % even spacing of the grid in x-dimension
end

% the grid is compressed in y-dimension, with depth Lmin at the
% shore and Lmax at the center of the lake.
for i=1:p.Zn
    for j = 1:p.Xn
        p.Z(i,j) = (p.Lmax/(p.Zn-1))*(i-1)*(1 + (p.Lmin/p.Lmax -1)* (p.X(i,j)/p.W).^(p.alpha));
    end
end

% coordinates of the center of each mesh quadrilateral
X_vol = zeros(p.Zn-1, p.Xn-1);
Z_vol = zeros(p.Zn-1, p.Xn-1);


for j=1:p.Xn-1
    for i =1:p.Zn-1
        X_vol(i,j) = (p.X(i,j) + p.X(i+1,j) + p.X(i,j+1) + p.X(i+1,j+1) ) /4;
        Z_vol(i,j) = (p.Z(i,j) + p.Z(i+1,j) + p.Z(i,j+1) + p.Z(i+1,j+1) ) /4;
    end
end
p.X_vol = X_vol;
p.Z_vol = Z_vol;

% [vol_areas, Area_bottom] = vol_areas_fn(p);
% p.Area_bottom = Area_bottom;
% p.vol_areas = vol_areas;

p.volumes_cyl = vol_areas_cyl_3d_fn(p);
p.Area_bottom_cyl = Area_bottom_cyl_fn(p);

%% Diffusion Coefficients
dx = zeros(p.Zn-1,p.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
dz = zeros(p.Zn-1,p.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]

diff_max_x = 1000; % horizontal diffusion coefficient at the surface
diff_max_z = 1000; % vertical diffusion coefficient at the surface
diff_min_x = 1000; % horizontal diffusion coefficient at the bottom
diff_min_z = 1000; % vertical diffusion coefficient at the bottom


% Linearly decreasing diffusion coefficients with the lake depth.
for j = 1:p.Xn-1
    for i = 1:p.Zn-1
        dx(i,j) =  diff_min_x + (diff_max_x-diff_min_x)*(1- p.Z_vol(i,j)/p.Lmax);
        dz(i,j) =  diff_min_z + (diff_max_z-diff_min_z)*(1- p.Z_vol(i,j)/p.Lmax);
    end
end

p.dx = dx;
p.dz = dz;
%  p.dx(1:end,:) = 1;
%  p.dx(1,:) = 10;
%  p.dx(2,:) = 9;
%  p.dx(3,:) = 8;
%  p.dx(4,:) = 7;
%  p.dx(5,:) = 6;
%  p.dx(6,:) = 5;
%  p.dx(7,:) = 4;
%  p.dx(8,:) = 3;
%  p.dx(9,:) = 2;
%  p.dx(10,:) = 1;

%% variable resuspension coefficient
% A type II functional response function is used to
%p.resus_depth

%% Pre-evaluation of transform derivatives
[dX_dXi_preCalc, dX_dEta_preCalc, dY_dXi_preCalc, dY_dEta_preCalc] = preCalc_Derivatives(p);
p.dX_dXi_preCalc  = dX_dXi_preCalc;
p.dX_dEta_preCalc  = dX_dEta_preCalc;
p.dZ_dXi_preCalc  = dY_dXi_preCalc;
p.dZ_dEta_preCalc = dY_dEta_preCalc;

%% Stiffness matrix & integration matrix
% the stiffness matrix is the matrix M in the linearized system y' = My,
% where y is a column vector.
p.S = Stiffness_matrix(p);
p.I_matrix = I_matrix_fn(p);

%% Inital Conditions

A0  = 1.0*ones(p.Zn-1, p.Xn-1);     % Initial Algal carbon density [mg C m^-3]
Rd0 = 10.000*ones(p.Zn-1, p.Xn-1);  % initial concentration of dissolved nutrients [mg P m^-3]
D0  = 1.000*ones(p.Zn-1, p.Xn-1);   % initial concentration of detritus [mg P m^-3]
Rs0 = 1.0*ones(1, p.Xn-1);          % initial concentration of sediment nutrient density [mg P m^-2]
B0  = 1.00*ones(1, p.Xn-1);         % initial concentration of benthic algal density [mg C m^-2]
  %A0(1,:)= 1.0;
%  A0(:,3)= 1.0;
 % A0(:,5)= 1.0;
%  A0(:,7)= 1.0;
  %A0(:,9)= 1.0;
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

 %Rd0(1,1) = 1;
% Rd0(1,2) = 1.1;
%Rd0(2,2) = 1.1;
% Rd0(2,2) = 2;
%Rd0(6,7) = 1;
%Rd0(16,16) = 5;

%Rd0(6,7) = 1;

%Rs0(1) = 1;

%B0(1) = 0;

%% calculation of total nutrient content at t=0, (to be conserved)
n_algae_0     = p.volumes_cyl.*A0*p.q;
n_dissolved_0 = p.volumes_cyl.*Rd0;
n_detritus_0  = p.volumes_cyl.*D0;
n_sediment_0  = Rs0.*p.Area_bottom_cyl;
n_benthic_0   = B0.*p.q_benth.*p.Area_bottom_cyl;

p.ntot_algae_0     = sum(sum(n_algae_0));
p.ntot_dissolved_0 = sum(sum(n_dissolved_0));
p.ntot_detritus_0  = sum(sum(n_detritus_0));
p.ntot_sed_0       = sum(n_sediment_0);
p.ntot_benthic_0   = sum(n_benthic_0);

ntot_0 =  p.ntot_algae_0 + p.ntot_dissolved_0 + p.ntot_detritus_0 + p.ntot_sed_0 + p.ntot_benthic_0;

%% state variables
A = A0'; % transpose of matrices in order to use colon notation to reshape to vector form.
Rd = Rd0';
D = D0';
Rs = Rs0';
B = B0';
y0 = [A(:); Rd(:); D(:); Rs(:); B(:)];

%% Simulation of model
%tend = 1.2518e+09; % end simulation time
% t_span = [1:tend/15: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.

%ode_opts = odeset( 'abstol' , 1e-7 , 'reltol', 1e-7, 'NonNegative',find(y0));
ode_opts = odeset( 'abstol' , 1e-9 , 'reltol', 1e-9, 'NonNegative',find(y0), 'Events', @(t,y) eventfun_V4(t,y,p));
%ode_opts = odeset(  'reltol', 1e-8, 'NonNegative',find(y0), 'Events', @(t,y) eventfun_V4(t,y,p));
[t,Y_t] = ode15s( @(t,Y) rhs_function(t,Y,p), [0, inf], y0, ode_opts);
%[t,Y_t] = ode15s( @(t,Y) dAdt_efficient_correct_V4_detritus(t,Y,p) ,[1 tend] , y0 , ode_opts);

%% recording of simulation time & saving workspace
simTime = toc;
%file_name = "2D_benthic_results_V4_dx_" + num2str(diff_max_z) + "_dz_" + num2str(diff_max_z) + "_Xn_" + num2str(p.Xn) + "_Zn_" + num2str(p.Zn);
%save(file_name);

%% Extracting results

A  = Y_t(end,1:(p.Xn-1)*(p.Zn-1));
Rd = Y_t(end,(p.Xn-1)*(p.Zn-1)+1 : 2*(p.Xn-1)*(p.Zn-1));
D  = Y_t(end,2*(p.Xn-1)*(p.Zn-1)+1 : 3*(p.Xn-1)*(p.Zn-1));
Rs = Y_t(end,3*(p.Xn-1)*(p.Zn-1) +1 : 3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1));
B  = Y_t(end,3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1)+1: end);
% Reshaping the vectors into matrices
A  = reshape(A, [p.Xn-1, p.Zn-1]);
Rd = reshape(Rd, [p.Xn-1, p.Zn-1]);
D  = reshape(D, [p.Xn-1, p.Zn-1]);
A  = A';
Rd = Rd';
D  = D';

%Calculating nutrient content at t = Tend
n_algae_end     = p.volumes_cyl.*A*p.q;
n_dissolved_end = p.volumes_cyl.*Rd;
n_detritus_end  = p.volumes_cyl.*D;
n_sediment_end  = Rs.*p.Area_bottom_cyl;
n_benthic_end   = B.*p.q_benth.*p.Area_bottom_cyl;
ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);

%% Calculation of the light intensity
I = zeros(p.Zn-1,p.Xn-1);

for x = 1:p.Xn-1
    integral = 0;
    for z=1:p.Zn-1
        integral = integral + p.Z(z,x) *(p.kA*A(z,x) + p.kD*D(z,x));
        I(z,x) = p.I0 * exp(-integral -p.kbg*p.Z(z,x));
    end
end

%% Figure index
figureIndex = 1;

%% Plot of simulation results
if(true)
    % creation of refined grid for interpolation.
    [X_vol_new,Y_vol_new] = grid_interpolation_fn(p);
    
   
    % plot of algal density
    if(true)
        % Plots algal density where the values have been interpolated to
        % produce a smooth figure.
        grid_data_A = griddata(p.X_vol,p.Z_vol, p.q.*A,X_vol_new,Y_vol_new); % interpolates to refined grid.
        figure(figureIndex)
        movegui(figureIndex,[1720 700]);
        figureIndex = figureIndex +1;
        h1 = axes;
        
        surf(X_vol_new,Y_vol_new,grid_data_A ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(h1, 'Ydir', 'reverse')
        %colorbar
        title("Concentration of Algae [mgC/m^3]");
        ylabel("Depth [m]");
        xlabel("distance from lake center [m]");
        zlabel("Algae concentration [mgC/m^3]");
        colorbar
        colormap(jet)
        %shading interp
    end
    
    % plot of dissolved nutrient density
    if(true)
        grid_data_Rd = griddata(p.X_vol,p.Z_vol, Rd,X_vol_new,Y_vol_new); % interpolates to refined grid.
        figure(figureIndex)
        movegui(figureIndex,[2300 700]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,grid_data_Rd ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colorbar
        colormap(jet)
        title("Concentration of dissolved nutrients [mgP/m^3]");
    end
    
    % plot of detritus nutrient density
    if(true)
        grid_data_D = griddata(p.X_vol,p.Z_vol, D, X_vol_new,Y_vol_new); % interpolates to refined grid.
        figure(figureIndex)
        movegui(figureIndex,[2880 700]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,grid_data_D ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colorbar
        colormap(jet)
        title("Concentration of detritus [mgP/m^3]");
    end
    
    % plot of sediment nutrient density
    if(true)
        % creating points for sediment plot
        sed_points = zeros(1,p.Xn-1);
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        
        figure(figureIndex)
        movegui(figureIndex,[1720 170]);
        figureIndex = figureIndex +1;
        x = [1:(p.W-1)/(p.Xn-1):p.W];
        plot(sed_points,Rs);
        title("Concentration of Sediment Nutrients [mgP/m^2]");
    end
    
    % plot of Benthic algal density [mgC/m^2]
    if(true)
        % creating points for benthic algae plot
        benth_points = zeros(1,p.Xn-1);
        benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        
        figure(figureIndex)
        movegui(figureIndex,[2300 170]);
        figureIndex = figureIndex +1;
        x = [1:(p.W-1)/(p.Xn-1):p.W];
        plot(sed_points,B);
        title("Concentration of Benthic Algae [mgC/m^2]");
    end
    
    %plot of algal production G
    if(false)
        figure(figureIndex)
        figureIndex = figureIndex +1;
        x0=400;
        y0=400;
        width=550;
        height=500;
        set(gcf,'position',[x0,y0,width,height])
        
        subaxis(1, 1, 1, 'sh', 0.000, 'sv', 0.00, ...
            'PL', 0.100,'PR', 0.10, 'PT', 0.10, 'PB', 0.10,...
            'MT', 0.0,'MB', 0.0, 'ML', 0.0, 'MR', 0.00);
        
        colormap(jet)
        %surf(p.X_vol,p.Z_vol ,griddata(X_volOld,Y_volOld, G,p.X_vol,p.Z_vol),  'edgecolor','none')
        surf(p.X_vol,p.Z_vol ,G,  'edgecolor','none')
        view(2);
        set(gca, 'YDir','reverse', 'XDir', 'reverse')
        axis square
        colorbar
        
        xticks(10:10:50);
        yticks(10:10:50);
        xlim([0 50])
        ylim([0 50])
        
        %%create a second axes.
        ax1 = gca; % the first axes
        ax1.FontSize = 28;
        ax1.FontWeight = 'bold';
        ax2 = axes('Position',ax1.Position,...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'Color','none',...
            'Ylim',ax1.YLim,...
            'XLim',ax1.XLim,...
            'TickLength',[0 0],...
            'YTick', [ax1.YLim(1):5:ax1.YLim(2)], ...
            'XTick', [ax1.XLim(1):5:ax1.XLim(2)],  ...
            'YTickLabel', [],  ...
            'XTickLabel', []  );
        ax2.FontSize = 28;
        linkaxes([ax1 ax2],'xy')
        grid on
        axis square
        temp = colorbar;
        set(temp,'YTick',[])
        set(gcf,'CurrentAxes',ax1); % return to ax1
        set(gcf,'CurrentAxes',ax1); % return to ax1
        cbar = colorbar;
        cbar.FontWeight = 'bold';
    end
    
    %plot of Algal light/nutrient limitation
    if(true)
        
        % calculation  of light/nutrient limitation at each point in lake.
        
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        
        
        interp_l_lim = griddata(p.X_vol,p.Z_vol, l_lim, X_vol_new,Y_vol_new); % interpolates to refined grid.
        interp_n_lim = griddata(p.X_vol,p.Z_vol, n_lim, X_vol_new,Y_vol_new); % interpolates to refined grid.
        
        limitation = zeros(length(interp_l_lim(:,1)),length(interp_l_lim(1,:)));
        for i = 1:length(interp_l_lim(:,1))
            for j =1:length(interp_l_lim(1,:))
                
                if(~isnan(interp_l_lim(i,j)) || ~isnan(interp_l_lim(i,j)))
                    if( interp_l_lim(i,j) > interp_n_lim(i,j))
                        limitation(i,j) = 1; % nutrient limited.
                    else
                        limitation(i,j) = -1; % light limited.
                    end
                else
                    limitation(i,j) = NaN;
                end
            end
        end
        
        
 
        figure(figureIndex)
        movegui(figureIndex,[2880 170]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,limitation ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
       % colorbar
        colormap(jet)
        title("light vs nutrient limitation");
        xlabel("red = nutrient limited, blue = light limited.");
    end
    
    % plot of lake coordinate system
    if(false)
        figure(figureIndex)
        figureIndex = figureIndex +1;
        h1 = axes;
        plot(p.X,p.Z,'*','Color', 'b');
        title('lake geometry');
        ylabel('depth [m]');
        xlabel('distance from lake center [m]');
        set(h1, 'Ydir', 'reverse')
        set(h1, 'YAxisLocation', 'Right')
    end
    
    % plot of initial concentration
    if(false)
        figure(figureIndex)
        figureIndex = figureIndex +1;
        h1 = axes;
        surf(p.X_vol,p.Z_vol,A0);
        title('initial Algae concentration');
        ylabel('depth [m]');
        xlabel('distance from lake center [m]');
        set(h1, 'Ydir', 'reverse')
        set(h1, 'YAxisLocation', 'Right')
    end
    
    
    perc_algae = sum(sum(n_algae_end))/ntot_end;
    perc_dissolved = sum(sum(n_dissolved_end))/ntot_end;
    perc_detritus = sum(sum(n_detritus_end))/ntot_end;
    perc_sediment = sum(n_sediment_end)/ntot_end;
    perc_benthic = sum(n_benthic_end)/ntot_end;
    
    disp("portion of leaked nutrients: " + num2str(abs((ntot_0-ntot_end) /ntot_0 )));
    disp("leaked nutrients [mgP]: " + num2str(abs(ntot_0-ntot_end)));
    disp("portion of nutrients in phytoplankton: " + num2str(perc_algae));
    disp("portion of nutrients in dissolved nutrients: " + num2str(perc_dissolved));
    disp("portion of nutrients in detritus: " + num2str(perc_detritus));
    disp("portion of nutrients in sediment: " + num2str(perc_sediment));
    disp("portion of nutrients in benthic algae: " + num2str(perc_benthic));
    disp("_____________________________________________");
end

%% plot of nutrient leakage over time
if(true)
leakage   = zeros(1,length(t));
algae     = zeros(1,length(t));
dissolved = zeros(1,length(t));
beth      = zeros(1,length(t));


for i = 1:length(t)
A  = Y_t(i,1:(p.Xn-1)*(p.Zn-1));
Rd = Y_t(i,(p.Xn-1)*(p.Zn-1)+1 : 2*(p.Xn-1)*(p.Zn-1));
D  = Y_t(i,2*(p.Xn-1)*(p.Zn-1)+1 : 3*(p.Xn-1)*(p.Zn-1));
Rs = Y_t(i,3*(p.Xn-1)*(p.Zn-1) +1 : 3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1));
B  = Y_t(i,3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1)+1: end);
% Reshaping the vectors into matrices
A  = reshape(A, [p.Xn-1, p.Zn-1]);
Rd = reshape(Rd, [p.Xn-1, p.Zn-1]);
D  = reshape(D, [p.Xn-1, p.Zn-1]);
A  = A';
Rd = Rd';
D  = D';

%Calculating nutrient content at t = Tend
n_algae_end     = p.volumes_cyl.*A*p.q;
n_dissolved_end = p.volumes_cyl.*Rd;
n_detritus_end  = p.volumes_cyl.*D;
n_sediment_end  = Rs.*p.Area_bottom_cyl;
n_benthic_end   = B.*p.q_benth.*p.Area_bottom_cyl;
ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);
leakage(i)      = (-ntot_end + ntot_0)/ntot_0;
end
figure(figureIndex)
figureIndex = figureIndex +1;
plot(t,leakage);
title('(ntot(t)-ntot0)/ ntot0');
end

%% video of plankton dynamics
A_vid = true; % flags for videos, set to true if desired.
Rd_vid = true;
Rs_vid = true;
Ab_vid = true;
D_vid = true;
lim_vid = true;
% creating points for sediment plot
sed_points = zeros(1,p.Xn-1);
sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;


if(false) % set to false if no videos are desired.
    close all
    if(A_vid) A_video = VideoWriter('plankton.avi');                open(A_video); end
    if(Rd_vid) Rd_video = VideoWriter('dissolved_nutrients.avi');   open(Rd_video);end
    if(D_vid) D_video = VideoWriter('dissolved_nutrients.avi');   open(D_video);end
    if(Rs_vid) Rs_video = VideoWriter('Sediment_nutrients.avi');    open(Rs_video); end
    if(Ab_vid) Ab_video = VideoWriter('benthic_algae.avi');         open(Ab_video); end
    if(lim_vid) lim_video = VideoWriter('benthic_algae.avi');         open(lim_video); end
    
    %axis tight manual
    set(gca,'nextplot','replacechildren');
    
    for t_index = 1:100:(length(Y_t(:,1)))
        % Extracting state variables
        if(A_vid)
            A = Y_t(t_index,1:(p.Xn-1)*(p.Zn-1));
            A = reshape(A, [(p.Xn-1), (p.Zn-1)]);
            A = A';
        end
        if(Rd_vid)
            Rd = Y_t(t_index, (p.Xn-1)*(p.Zn-1)+1 : 2*(p.Xn-1)*(p.Zn-1));
            Rd = reshape(Rd, [(p.Xn-1), (p.Zn-1)]);
            Rd = Rd';
        end
        if(D_vid)
            D = Y_t(t_index, 2*(p.Xn-1)*(p.Zn-1)+1 : 3*(p.Xn-1)*(p.Zn-1));
            D = reshape(D, [(p.Xn-1), (p.Zn-1)]);
            D = D';
        end
        if(Rs_vid)
            Rs = Y_t(t_index, 3*(p.Xn-1)*(p.Zn-1) +1 : 3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1));
        end
        if(Ab_vid)
            B = Y_t(t_index, 3*(p.Xn-1)*(p.Zn-1) + (p.Xn-1) +1 : end);
        end
        
        colormap(jet)
        %%%% Nutrients bound in Algae %%%%
        if(A_vid)
            f1 = figure(1);
            movegui(f1,[1720 700]);
            clf(f1);
            
            
            h6 = axes;
            %grid_data = griddata(p.X_vol,p.Z_vol, A,X_vol_new,Y_vol_new);
      % 
            %surf(X_vol_new,Y_vol_new,grid_data ,  'edgecolor','none')
            %caxis([0,max(max(A0))]);
            surf(p.X_vol, p.Z_vol, A, 'edgecolor','none');
            shading interp
            grid off
            ylabel("Depth [m]");
            xlabel("distance from lake center [m]");
            title("Algae concentration [mgC/m^3]");
            colorbar;
            az =0;
            el = 90;
            %az= 45+180; % azimuth angle
            %el= 65; % elevation angle
            view(az,el);
            set(h6, 'Ydir', 'reverse')
            %axis([0 p.W 0 p.Lmax 0 max(max(A))]) % Fixing window size for video
            
            ax = struct('Axes', gca);
            %align_axislabel([],ax) % aligns the axis labels to the plot
            
            if(az ~=0 && el ~= 90)
                align_axislabel([],ax) % aligns the axis labels to the plot
            end
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(A_video,frame);
        end
        
        %%%% Dissolved Nutrients in the Water %%%%
        if(Rd_vid)
            f2 = figure(2);
            movegui(f2,[2300 700]);
            clf(f2);
            h6 = axes;
            %grid_data = griddata(p.X_vol,p.Z_vol, Rd,X_vol_new,Y_vol_new);
            %surf(X_vol_new,Y_vol_new,grid_data ,  'edgecolor','none')
            surf(p.X_vol, p.Z_vol, Rd,'edgecolor','none');
            shading interp
            grid off
            ylabel("Depth [m]");
            xlabel("distance from lake center [m]");
            title("Dissolved Nutrient concentration [mgP/m^3]");
            colorbar;
            az =0;
            el = 90;
            %az= 45+180; % azimuth angle
            %el= 65; % elevation angle
            view(az,el);
            set(h6, 'Ydir', 'reverse')
            %axis([0 p.W 0 p.Lmax 0 max(max(A0))]) % Fixing window size for video
            %axis([0 p.W 0 p.Lmax 0 2]) % Fixing window size for video
            %ax = struct('Axes', gca);
            %align_axislabel([],ax) % aligns the axis labels to the plot
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(Rd_video,frame);
        end
        
        %%%% Detritus in the water %%%%
        if(D_vid)
            f3 = figure(3);
            movegui(f3,[2880 700]);
            clf(f3);
            h6 = axes;
            %grid_data = griddata(p.X_vol,p.Z_vol, D,X_vol_new,Y_vol_new);
            %surf(X_vol_new,Y_vol_new,grid_data ,  'edgecolor','none')
            surf(p.X_vol, p.Z_vol, D,'edgecolor','none');
            shading interp
            grid off
            ylabel("Depth [m]");
            xlabel("distance from lake center [m]");
            title("Detritus Nutrient concentration [mgP/m^3]");
            colorbar;
            az =0;
            el = 90;
            %az= 45+180; % azimuth angle
            %el= 65; % elevation angle
            view(az,el);
            set(h6, 'Ydir', 'reverse')
            %axis([0 p.W 0 p.Lmax 0 max(max(A0))]) % Fixing window size for video
            %axis([0 p.W 0 p.Lmax 0 2]) % Fixing window size for video
            %ax = struct('Axes', gca);
            %align_axislabel([],ax) % aligns the axis labels to the plot
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(D_video,frame);
        end
        
        %%%% Nutrients in the Sediment %%%%
        if(Rs_vid)
            f4 = figure(4);
            movegui(f4,[1720 170]);
            clf(f4);
            %h6 = axes;
            Rs_interp = interp1(sed_points,Rs, linspace(0,p.W,300));
            plot(linspace(0,p.W,300),Rs_interp);
            %plot(sed_points,Rs);
            grid off
            ylabel("Sediment nutrient concentration [mgP/m^2]");
            xlabel("distance from lake center [m]");
            title("Sediment nutrient concentration [mgP/m^2]");
            %set(h6, 'Ydir', 'reverse')
            %axis([0 p.W 0 p.Lmax 0 max(max(A0))]) % Fixing window size for video
            %axis([0 p.W 0 p.Lmax -2 2]) % Fixing window size for video
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(Rs_video,frame);
        end
        
        %%%% Benthic algae [mgC/m^2] %%%%
        if(Ab_vid)
            f5 = figure(5);
            movegui(f5,[2300 170]);
            clf(f5);
            %plot(sed_points,B);
            B_interp = interp1(sed_points,B, linspace(0,p.W,300));
            plot(linspace(0,p.W,300),B_interp);
            grid off
            ylabel("Benthic Algae concentration [mgC/m^2]");
            xlabel("distance from lake center [m]");
            title("Benthic Algae concentration [mgC/m^2]");
            drawnow
            frame = getframe(gcf);
            writeVideo(Ab_video,frame);
        end
        
        %%%% light/nutrient limitation plot %%%%
        if(lim_vid)
             % calculation  of light/nutrient limitation at each point in lake.
        
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        
        
        interp_l_lim = griddata(p.X_vol,p.Z_vol, l_lim, X_vol_new,Y_vol_new); % interpolates to refined grid.
        interp_n_lim = griddata(p.X_vol,p.Z_vol, n_lim, X_vol_new,Y_vol_new); % interpolates to refined grid.
        
        limitation = zeros(length(interp_l_lim(:,1)),length(interp_l_lim(1,:)));
        for i = 1:length(interp_l_lim(:,1))
            for j =1:length(interp_l_lim(1,:))
                
                if(~isnan(interp_l_lim(i,j)) || ~isnan(interp_l_lim(i,j)))
                    if( interp_l_lim(i,j) > interp_n_lim(i,j))
                        limitation(i,j) = 1; % nutrient limited.
                    else
                        limitation(i,j) = -1; % light limited.
                    end
                else
                    limitation(i,j) = NaN;
                end
            end
        end
        
        f6 = figure(6);
        movegui(f6,[2880 170]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,limitation ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colorbar
        colormap(jet)
        title("light vs nutrient limitation");
        xlabel("red = nutrient limited (1), blue = light limited (-1).");
        end
        
    end
    
    
    if(A_vid)  close(A_video);  end
    if(Rd_vid) close(Rd_video); end
    if(Rs_vid) close(Rs_video); end
    if(Ab_vid) close(Ab_video); end
    if(lim_vid) close(lim_video); end
end