
%% Implementation details and a clean slate
close all
clear
clc
tic
% This code implements the system of equations described in Jäger et al. (2010)
% in one dimension (z-dim) using a finite difference method, with an added
% layer of benthic algae residing on top of the sediment layer.
% The lake is assumed to be a cylindrical shell with an inner radius R_inner, outer
% radius R_outer, and constant depth Lmax when calculating the nutrient
% content. This is done because the implementation is intended to be used in
% comparison with a 2d cylindrical implementation of the same model.
% ode15s is used to iterate in time. Plankton nutrient density stochiometry
% q is constant.
% 
% Written by Hugo Harlin 2019

%% Model Parameters
% all parameters are organized in a struct p for convenience
p = struct;
p.dz = 10.0; % Radial Turbulent-diffusion coefficient [m^2 day^-1]
p.I0 = 300; % Light intensity at the surface [micro-mol photons m^-2 s^-1]
p.k = 0.0003; % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
p.kB = 0.0003; % Specific light-attenuation coefficient of Benthic biomass [m^2 mg C^-1]
p.kbg = 0.4; % Background light-attenuation coefficient [m^-1]
p.lbg = 0.1; % Specific algal maintenance respiration losses [day^-1]
p.M = 1.5; % Half-saturation constant of algal nutrient uptake [mg P m^-3]
p.M_benth = 120.0; % Half-saturation constant of benthic algae nutrient uptake [mgPm^-3]
p.Gmax = 1.08; % Maximum specific phytoplankton production rate [day^-1]
p.Gmax_benth = 1.08; % Maximum specific benthic algae production rate [day^-1]
p.H = 120.0; % Half-saturation constant of light-dependent algal production [micro-mol photons m^-2 s^-1]
p.q = 0.00244; % Algal nutrient quota, Redfield ratio [mgP/mgC]
p.q_benth = 0.0244; % Benthic algae nutrient quota, Redfield ratio [mgP/mgC]
p.lbg_benth = 0.01; % Specific benthic algae maintenance respiration losses [day^-1]
p.r = 0.02; % Specific remineralization rate of sedimented nutrients [day^-1]
p.v = 0.0; % Algal sinking velocity [m day^-1]
%Not implemented: p.death_rate = 1; % coefficient governing the proportion of sinking algae at the bottom that dies.
                  % 0 = no death. 1 = all algae that would have sunk through the sediment dies. 
p.benth_recycling = 1.0; % Governs the portion of respired nutrients that are released as dissolved nutrients.
                         % the rest is bound in partculate matter in the sediment. 
                       
%% Lake topology and Mesh
% Quantities relating to system size. The script is written to facilitate
% iteration over the radii and the depth.
p.X_inner = 0; % inner radius
p.X_outer = 1; % outer radius
p.Zn = 17; % Number of grid-points (width)
p.Lmax = 10; % lake depth [m]
p.deltaZ = p.Lmax/(p.Zn-1);
p.Z_nodes = linspace(0,p.Lmax, p.Zn); % Mesh nodes

% Cell centered mesh
for i=1:p.Zn-1
    p.Z(i) = (p.Z_nodes(i)+p.Z_nodes(i+1))/2;
end

%% Inital Conditions

A0 = 0.001*ones(1,p.Zn-1); % Initial Algal carbon density [mg C m^-3]
A0(5) = 0;

Rd0 = 1.00*ones(1,p.Zn-1); % initial concentration of dissolved nutrients [mg P m^-3]
Rd0(8) = 0;

Rs0 = 1; % initial concentration of sediment nutrient density [mg P m^-2]
%Rs0 = 1;
B0 = 1; % initial concentration of benthic algal density [mg C m^-2]

% calculation of total nutrient content at t=0, to be conserved
n_algae_0     = p.Lmax/(p.Zn-1)*pi*(p.X_outer^2-p.X_inner^2)*A0*p.q;
n_dissolved_0 = p.Lmax/(p.Zn-1)*pi*(p.X_outer^2-p.X_inner^2)*Rd0;
n_sediment_0  = Rs0.*pi*(p.X_outer^2-p.X_inner^2);
n_benthic_0   = B0.*p.q_benth*pi*(p.X_outer^2-p.X_inner^2);

p.ntot_algae_0  = sum(n_algae_0);
p.ntot_dissolved_0 = sum(n_dissolved_0);
p.ntot_sed_0 = n_sediment_0;
p.ntot_benthic_0 = n_benthic_0;

ntot_0 =  p.ntot_algae_0 + p.ntot_dissolved_0 + p.ntot_sed_0 + p.ntot_benthic_0;

%% state variables
A = A0; % transpose of matrices in order to use colon notation to reshape to vector form.
Rd = Rd0;
Rs = Rs0;
B = B0;
y0 = [A(:); Rd(:); Rs(:); B(:)];
I = zeros(1,p.Zn);

%% matrices for finite difference method

p.diffusion_matrix =  diag(-2*ones(1,p.Zn-1)) + diag(1*ones(1,p.Zn-2),1) + diag(1*ones(1,p.Zn-2),-1);
p.diffusion_matrix(1,1) = -1;
p.diffusion_matrix(end,end) = -1;
p.diffusion_matrix = p.diffusion_matrix.*p.dz/(p.deltaZ^2);


p.sink_matrix =  diag(-1*ones(1,p.Zn-1)) + diag(1*ones(1,p.Zn-2),-1);
p.sink_matrix = p.sink_matrix*p.v/p.deltaZ;
%p.sink_matrix(1,1) = -1*p.sink_matrix(1,1);
%p.sink_matrix(end,end) = 0;

p.dynam_matrix = p.sink_matrix + p.diffusion_matrix;
%% Simulation of model
tend = 400; % end simulation time
t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
%t_span = [1:tend];
miac = @(t,y) eventfun(t,y,p);
ode_opts = odeset('NonNegative',[1:2*(p.Zn-1)+2] , 'abstol' , 1e-8 , 'reltol' , 1e-8,'Events',miac);
%ode_opts = odeset( 'abstol' , 1e-6 , 'reltol' , 1e-9);
[t,Y_t] = ode15s( @(t,Y) dAdt(t,Y,p,I) , [0 inf] , y0 , ode_opts); 
%[t,Y_t] = ode15s( @(t,Y) dAdt(t,Y,p,I) , t_span , y0 , ode_opts);

%% recording of simulation time & saving workspace
simTime = toc;
file_name = "_dz_" + num2str(p.dz) + "1D" ;
save(file_name);

%% Extracting results

A = Y_t(end, 1:(p.Zn-1));
Rd = Y_t(end, p.Zn:2*(p.Zn-1));
Rs = Y_t(end, 2*(p.Zn-1)+1);
B = Y_t(end, 2*(p.Zn-1)+2);

%Calculating nutrient content at t = Tend
n_algae_end = pi*p.Lmax/(p.Zn-1)*(p.X_outer^2-p.X_inner^2)*A*p.q;
n_dissolved_end = pi*p.Lmax/(p.Zn-1)*(p.X_outer^2-p.X_inner^2)*Rd;
n_sediment_end = Rs.*pi*(p.X_outer^2-p.X_inner^2);
n_benthic_end = B.*p.q_benth*pi*(p.X_outer^2-p.X_inner^2);
ntot_end =  sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(n_sediment_end) + sum(n_benthic_end);

% Calculates the algal production G at each point in the computation grid
% and returns it as a matix

% Calculation of the light intensity at each grid element center
% I = zeros(1,p.Zn);

% Calculation of the algal production G(R,I)
%G =  p.Gmax .* Rd./(Rd+p.M) .* I./(I+p.H);
%GA = G.*A;

%% Figure index
figureIndex = 1;

%% Plot of simulation results
close all
figureIndex = 1;
colormap(jet)
% plot of algal density
if(true)
    % Plots algal density where the values have been interpolated to
    % produce a smooth figure.
    figure(figureIndex)
    movegui(figureIndex,[1750 700]);
    figureIndex = figureIndex +1;
    h1 = axes;
    
    plot(p.Z,A);
    title("Concentration of Algae [mgC/m^3]");
    ylabel("Algae concentration [mgC/m^3]");
    xlabel("depth [m]");
    colorbar
    %shading interp
end

% plot of dissolved nutrient density
if(true)
    figure(figureIndex)
    movegui(figureIndex,[2350 700]);
    figureIndex = figureIndex +1;
    plot(p.Z,Rd);
    colorbar
    colormap(jet)
    title("Concentration dissolved nutrients [mgP/m^3]");
end

% plot of sediment nutrient density
if(false)
    % creating points for sediment plot
    sed_points = zeros(1,p.Xn-1);
    sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
    
    figure(figureIndex)
    movegui(figureIndex,[1750 170]);
    figureIndex = figureIndex +1;
    x = [1:(p.W-1)/(p.Xn-1):p.W];
    plot(sed_points,Rs);
    title("Concentration of Sediment Nutrients [mgP/m^2]");
end

% plot of Benthic algal density [mgC/m^2]
if(false)
    % creating points for benthic algae plot
    benth_points = zeros(1,p.Xn-1);
    benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
    
    figure(figureIndex)
    movegui(figureIndex,[2350 170]);
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

%plot of net algal production G*A
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
    %surf(p.X_vol,p.Z_vol ,griddata(X_volOld,Y_volOld, GA,p.X_vol,p.Z_vol),  'edgecolor','none')
    surf(p.X_vol,p.Z_vol,GA,  'edgecolor','none')
    view(2);
    set(gca, 'YDir','reverse', 'XDir', 'reverse')
    axis square
    %colorbar
    
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

%Light intensity plot at the points (x_coord,:)
if(false)
    figure(figureIndex)
    figureIndex = figureIndex +1;
    
    x_coord = 1;
    a = zeros(1,p.Zn);
    depth = p.Z(:,x_coord);
    for i=1:p.Zn
        a(i) = light_intensity(p,A,x_coord,i);
    end
    plot(depth,a)
    title('Light intensity' );
    ylabel('Light Intensity [micro-mol photons m^-1 s^-1]');
    xlabel('Water depth [m]');
end

perc_algae = sum(sum(n_algae_end))/ntot_end;
perc_dissolved = sum(sum(n_dissolved_end))/ntot_end;
perc_sediment = sum(n_sediment_end)/ntot_end;
perc_benthic = sum(n_benthic_end)/ntot_end;

disp("percentage leaked nutrients: " + num2str(abs((ntot_0-ntot_end) /ntot_0 )));
disp("portion of nutrients in phytoplankton: " + num2str(perc_algae));
disp("portion of nutrients in dissolved nutrients: " + num2str(perc_dissolved));
disp("portion of nutrients in sediment: " + num2str(perc_sediment));
disp("portion of nutrients in benthic algae: " + num2str(perc_benthic));
disp("_____________________________________________");


%% video of plankton dynamics
A_vid = false; % flags for videos, set to true if desired.
Rd_vid = false;
Rs_vid = false;
Ab_vid = false;
% creating points for sediment plot
%sed_points = zeros(1,p.Xn-1);
%sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;


if(false) % set to false if no videos are desired.
    close all
    if(A_vid) A_video = VideoWriter('plankton.avi');                open(A_video); end
    if(Rd_vid) Rd_video = VideoWriter('dissolved_nutrients.avi');   open(Rd_video);end
    if(Rs_vid) Rs_video = VideoWriter('Sediment_nutrients.avi');    open(Rs_video); end
    if(Ab_vid) Ab_video = VideoWriter('benthic_algae.avi');         open(Ab_video); end
    %axis tight manual
    set(gca,'nextplot','replacechildren');
    
    for t = 1:1:length(Y_t(:,1))
        % Extracting state variables
        if(A_vid)
            A = Y_t(t,1:(p.Zn-1));
        end
        if(Rd_vid)
            Rd = Y_t(t, (p.Zn-1)+1 : 2*(p.Zn-1));
        end
        if(Rs_vid)
            Rs = Y_t(t, 2*(p.Zn-1) +1 : 2*(p.Zn-1) + (p.Zn-1));
        end
        if(Ab_vid)
            B = Y_t(t, 2*(p.Xn-1)*(p.Zn-1) + (p.Xn-1) +1 : end);
        end
        
        colormap(jet)
        %%%% Nutrients bound in Algae %%%%
        if(A_vid)
            f1 = figure(1);
            movegui(f1,[1750 700]);
            clf(f1);
            
            plot(p.Z,A)
            grid off
            ylabel("Algae concentration [mgC/m^3]");
            xlabel("Depth [m]");
            title("Algae concentration [mgC/m^3]");
            colorbar;
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(A_video,frame);
        end
        
        %%%% Dissolved Nutrients in the Water %%%%
        if(Rd_vid)
            f2 = figure(2);
            movegui(f2,[2350 700]);
            clf(f2);
            plot(p.Z,Rd)
            grid off
            ylabel("Dissolved nutrient concentration [mgC/m^3]");
            xlabel("Depth [m]");
            title("Dissolved nutrient concentration [mgC/m^3]");
            colorbar;
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(Rd_video,frame);
        end
        
        %%%% Nutrients in the Sediment %%%%
        if(Rs_vid)
            f3 = figure(3);
            movegui(f3,[1750 170]);
            clf(f3);
            %h6 = axes;
            plot(sed_points,Rs);
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
            f4 = figure(4);
            movegui(f4,[2350 170]);
            clf(f4);
            %h6 = axes;
            plot(sed_points,B);
            grid off
            ylabel("Benthic Algae concentration [mgC/m^2]");
            xlabel("distance from lake center [m]");
            title("Benthic Algae concentration [mgC/m^2]");
            %set(h6, 'Ydir', 'reverse')
            %axis([0 p.W 0 p.Lmax 0 max(max(A0))]) % Fixing window size for video
            %axis([0 p.W 0 p.Lmax -2 2]) % Fixing window size for video
            drawnow
            %pause(0.0001)
            frame = getframe(gcf);
            writeVideo(Ab_video,frame);
        end        
        
        
    end
    if(A_vid)  close(A_video);  end
    if(Rd_vid) close(Rd_video); end
    if(Rs_vid) close(Rs_video); end
    if(Ab_vid) close(Ab_video); end    
end
