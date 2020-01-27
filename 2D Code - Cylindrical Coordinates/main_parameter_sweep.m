
alpha = 0.5;
for dx=[1,10,100,1000]
    for dy = [1,10,100,1000]
        dx
        dy
        %% Implementation details and a clean slate
        tic
        % This code implements the system of equations described in Jäger et al. (2010)
        % in two spatial dimensions, using finite volumes and Gauss divergence
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
        % Written by Hugo Harlin 2019
        
        %% Model Parameters
        % all parameters are organized in a struct p for implementation convenience
        p = struct;
        p.dx = dx; % Radial Turbulent-diffusion coefficient [m^2 day^-1]
        p.dy = dy; % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
        p.I0 = 300; % Light intensity at the surface [micro-mol photons m^-2 s^-1]
        p.k = 0.0003; % Specific light-attenuation coefficient of algal biomass [m^2 mg C^-1]
        p.kbg = 0.4; % Background light-attenuation coefficient [m^-1]
        p.lbg = 0.1; % Specific algal maintenance respiration losses [day^-1]
        p.M = 1.5; % Half-saturation constant of algal nutrient uptake [mg P m^-3]
        p.Gmax = 1.08; % Maximum specific algal production rate in the fixed-stochiometry and light-limitation models [day^-1]
        p.H = 120.0; % Half-saturation constant of light-dependent algal production [micro-mol photons m^-2 s^-1]
        p.q = 0.0244; % Algal nutrient quota, Redfield ratio
        p.r = 0.02; % Specific mineralization rate of sedimented nutrients [day^-1]
        p.v = 0.25; % Algal sinking velocity [m day^-1]
        
        %% Lake topology and Mesh
        
        cyl_cord = true;
        
        % Quantities relating to system size
        p.Xn = 30; % Number of grid-points (width)
        p.Yn = 30; % Number of grid-points (depth)
        p.Lmin = 0.00001; % Minimum lake depth (depth at land-water interface) [m]
        p.Lmax = 50; % Maximum lake depth [m]
        p.W = 50; % Lake radius [m]
        p.alpha = alpha; % Exponent governing the slope of the lake bottom
        
        % Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
        % at the center of the lake (slope = alpha* (Lmin - Lmax)/W).
        % (0,0) is placed at the center of the lake at the surface, y-dim is facing
        % downward and x-dim is facing towards the lake edge
        
        p.X = zeros(p.Yn, p.Xn); % Mesh-spacing in x-dimension
        p.Y = zeros(p.Yn, p.Xn); % Mesh-spacing in y-dimension
        
        for i=1:1:p.Yn
            p.X(i,:) = [0:p.W/(p.Xn-1):p.W]; % even spacing of the grid in x-dimension
        end
        
        % the grid is compressed in y-dimension, with depth Lmin at the
        % shore and Lmax at the center of the lake.
        for i=1:p.Yn
            for j = 1:p.Xn
                p.Y(i,j) = (p.Lmax/(p.Yn-1))*(i-1)*(1 + (p.Lmin/p.Lmax -1)* (p.X(i,j)/p.W).^(p.alpha));
            end
        end
        
        % coordinates of the center of each mesh quadrilateral
        [X_vol,Y_vol] = vol_element_coords(p);
        p.X_vol = X_vol;
        p.Y_vol = Y_vol;
        
        %% Pre-evaluation o f transform derivatives
        [dX_dXi_preCalc, dX_dEta_preCalc, dY_dXi_preCalc, dY_dEta_preCalc] = preCalc_Derivatives(p);
        p.dX_dXi_preCalc  = dX_dXi_preCalc;
        p.dX_dEta_preCalc  = dX_dEta_preCalc;
        p.dY_dXi_preCalc  = dY_dXi_preCalc;
        p.dY_dEta_preCalc = dY_dEta_preCalc;
        
        %% Inital Conditions
        
        % [vol_areas, L_bottom] = vol_areas_fn(p);
        % p.L_bottom = L_bottom;
        % p.vol_areas = vol_areas;
        
        p.volumes_cyl = vol_areas_cyl_3d_fn(p);
        p.L_bottom_cyl = L_bottom_cyl_fn(p);
        
        A0 = 1.00*ones(p.Yn-1, p.Xn-1); % Initial Algal carbon density [mg C m^-3]
        %A0(1,10)= 1;
        %A0(1,1)= 30;
        %A0(1,:)= 1;
        %A0(1,:) = 10./p.volumes_cyl(1,:);
        %A0(1,:) = 10;
        %A0(15,15)=10;
        %A0(:,p.Xn-4) = 10;
        Rd0 = 1.00*ones(p.Yn-1, p.Xn-1); % initial concentration of dissolved nutrients [mg P m^-3]
        %Rd0(end,:) = 0;
        Rs0 = 1*ones(1, p.Xn-1); % initial concentration of sediment nutrient density [mg P m^-2]
        %Rs0(1) = 1;
        
        if(~cyl_cord) % old function, cartesian coordinates.
            % calculation of total nutrient content at t=0, to be conserved
            [n_algae_0,n_dissolved_0,n_sediment_0, vol_areas, L_bottom] = Nutrient_content(A0,Rd0,Rs0,p);
            
            p.L_bottom = L_bottom;
            p.vol_areas = vol_areas;
            p.ntot_algae_0  = sum(sum(n_algae_0));
            p.ntot_dissolved_0 = sum(sum(n_dissolved_0));
            p.ntot_sed_0 = sum(n_sediment_0);
            
            ntot_0 =  p.ntot_algae_0 + p.ntot_dissolved_0 + p.ntot_sed_0;
        end
        
        if(cyl_cord)
            % calculation of total nutrient content at t=0, to be conserved
            n_algae_0 = p.volumes_cyl.*A0*p.q;
            n_dissolved_0 = p.volumes_cyl.*Rd0;
            n_sediment_0 = Rs0.*p.L_bottom_cyl;
            %
            %p.L_bottom = L_bottom;
            p.ntot_algae_0  = sum(sum(n_algae_0));
            p.ntot_dissolved_0 = sum(sum(n_dissolved_0));
            p.ntot_sed_0 = sum(n_sediment_0);
            
            ntot_0 =  p.ntot_algae_0 + p.ntot_dissolved_0 + p.ntot_sed_0;
        end
        
        %% state variables
        A = A0'; % transpose of matrices in order to use colon notation to reshape to vector form.
        Rd = Rd0';
        Rs = Rs0';
        y0 = [A(:); Rd(:); Rs(:)];
        I = zeros(p.Yn, p.Xn);
        
        %% Simulation of model
        tend = 5000; % end simulation time
        t_span = [1:tend/30: tend]; % timespan of simulation. The intermediate steps tells ode15s when to save current state of the model.
        %t_span = [1:tend];
        miac = @(t,y) eventfun(t,y,p);
        ode_opts = odeset( 'abstol' , 1e-6 , 'reltol' , 1e-9,'Events',miac);
        %ode_opts = odeset( 'abstol' , 1e-6 , 'reltol' , 1e-9);
        %[t,Y_t] = ode15s( @(t,Y) dAdt_fun(t,Y,p,I) , t_span , y0 , ode_opts ); %bug somewhere.
        [t,Y_t] = ode15s( @(t,Y) dAdt_alt(t,Y,p,I) , t_span , y0 , ode_opts ); % diffusion works in this version!
        
        %% recording of simulation time & saving workspace
        simTime = toc;
        %file_name = "dx_" + num2str(p.dx) + "_dz_" + num2str(p.dy) + "_alpha_" + num2str(p.alpha) ;
        file_name = "dx_" + num2str(p.dx) + "_dz_" + num2str(p.dy) + "_alpha_" + "0_5" ;
        save(file_name);
        
        %% Extracting results
        A = Y_t(end,1:(p.Xn-1)*(p.Yn-1));
        Rd = Y_t(end, (p.Xn-1)*(p.Yn-1)+1 : 2*(p.Xn-1)*(p.Yn-1));
        Rs = Y_t(end, 2*(p.Xn-1)*(p.Yn-1) +1 : 2*(p.Xn-1)*(p.Yn-1) + (p.Xn-1));
        
        % Reshaping the vectors into matrices
        A = reshape(A, [p.Xn-1, p.Yn-1]);
        Rd = reshape(Rd, [p.Xn-1, p.Yn-1]);
        A = A';
        Rd = Rd';
        
        if(~cyl_cord) % Incorrect for cylindrical coordiantes
            %Calculating nutrient content at t = Tend
            [n_algae_end,n_dissolved_end, n_sediment_end, vol_areas, L_bottom]  = Nutrient_content(A,Rd,Rs,p);
            ntot_end =  sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(n_sediment_end);
        end
        
        if(cyl_cord) % cylindrical coordinates
            %Calculating nutrient content at t = Tend
            n_algae_end = p.volumes_cyl.*A*p.q;
            n_dissolved_end = p.volumes_cyl.*Rd;
            n_sediment_end = Rs.*p.L_bottom_cyl;
            ntot_end =  sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(n_sediment_end);
        end
        
        % Calculating production G and gross production G*A
        G = G_fun(p,A,Rd);
        GA = G.*A;
        
        %% Figure index
        figureIndex = 1;
        
        %% Plot of simulation results
        close all
        figureIndex = 1;
        
        [X_vol_new,Y_vol_new] = grid_interpolation_fn(p,A);
        grid_data_A = griddata(p.X_vol,p.Y_vol, p.q.*A,X_vol_new,Y_vol_new);
        grid_data_Rd = griddata(p.X_vol,p.Y_vol, p.q.*Rd,X_vol_new,Y_vol_new);
        
        colormap(jet)
        % plot of algal density
        if(true)
            % Plots algal density where the values have been interpolated to
            % produce a smooth figure.
            figure(figureIndex)
            movegui(figureIndex,[1750 700]);
            figureIndex = figureIndex +1;
            h1 = axes;
            
            surf(X_vol_new,Y_vol_new,grid_data_A ,  'edgecolor','none')
            
            grid off
            az =0;
            el = 90;
            view(az,el);
            set(h1, 'Ydir', 'reverse')
            %colorbar
            title("Concentration of Algae");
            ylabel("Depth [m]");
            xlabel("distance from lake center [m]");
            zlabel("Algae concentration [mgC/m^3]");
            colorbar
            %shading interp
        end
        
        % plot of dissolved nutrient density
        if(true)
            figure(figureIndex)
            movegui(figureIndex,[2350 700]);
            figureIndex = figureIndex +1;
            surf(X_vol_new,Y_vol_new,grid_data_Rd ,  'edgecolor','none')
            grid off
            az =0;
            el = 90;
            view(az,el);
            set(gca, 'Ydir', 'reverse')
            colorbar
            colormap(jet)
            title("Concentration of dissolved nutrients");
        end
        
        % plot of sediment nutrient density
        if(true)
                    % creating points for sediment plot
        sed_points = zeros(1,p.Xn-1);
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        
            figure(figureIndex)
            movegui(figureIndex,[1750 170]);
            figureIndex = figureIndex +1;
            x = [1:(p.W-1)/(p.Xn-1):p.W];
            plot(sed_points,Rs);
            title("Concentration of Sediment Nutrients");
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
            %surf(p.X_vol,p.Y_vol ,griddata(X_volOld,Y_volOld, G,p.X_vol,p.Y_vol),  'edgecolor','none')
            surf(p.X_vol,p.Y_vol ,G,  'edgecolor','none')
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
            %surf(p.X_vol,p.Y_vol ,griddata(X_volOld,Y_volOld, GA,p.X_vol,p.Y_vol),  'edgecolor','none')
            surf(p.X_vol,p.Y_vol,GA,  'edgecolor','none')
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
            plot(p.X,p.Y,'*','Color', 'b');
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
            surf(p.X_vol,p.Y_vol,A0);
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
            a = zeros(1,p.Yn);
            depth = p.Y(:,x_coord);
            for i=1:p.Yn
                a(i) = light_intensity(p,A,x_coord,i);
            end
            plot(depth,a)
            title('Light intensity' );
            ylabel('Light Intensity [micro-mol photons m^-1 s^-1]');
            xlabel('Water depth [m]');
        end
        
        %% video of plankton dynamics
        A_vid = true; % flags for videos, set to true if desired.
        Rd_vid = true;
        Rs_vid = true;
        % creating points for sediment plot
        sed_points = zeros(1,p.Xn-1);
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        
        
        if(true) % set to false if no videos are desired.
            close all
            if(A_vid) A_video = VideoWriter('plankton.avi');                open(A_video); end
            if(Rd_vid) Rd_video = VideoWriter('dissolved_nutrients.avi');   open(Rd_video);end
            if(Rs_vid) Rs_video = VideoWriter('Sediment_nutrients.avi');    open(Rs_video); end
            
            %axis tight manual
            set(gca,'nextplot','replacechildren');
            
            for t = 1:1:length(Y_t(:,1))
                % Extracting state variables
                if(A_vid)
                    A = Y_t(t,1:(p.Xn-1)*(p.Yn-1));
                    A = reshape(A, [(p.Xn-1), (p.Yn-1)]);
                    A = A';
                end
                if(Rd_vid)
                    Rd = Y_t(t, (p.Xn-1)*(p.Yn-1)+1 : 2*(p.Xn-1)*(p.Yn-1));
                    Rd = reshape(Rd, [(p.Xn-1), (p.Yn-1)]);
                    Rd = Rd';
                end
                if(Rs_vid)
                    Rs = Y_t(t, 2*(p.Xn-1)*(p.Yn-1) +1 : 2*(p.Xn-1)*(p.Yn-1) + (p.Xn-1));
                end
                
                colormap(jet)
                %%%% Nutrients bound in Algae %%%%
                if(A_vid)
                    f1 = figure(1);
                    movegui(f1,[1750 700]);
                    clf(f1);
                    
                    
                    h6 = axes;
                    grid_data = griddata(p.X_vol,p.Y_vol, A,X_vol_new,Y_vol_new);
                    surf(X_vol_new,Y_vol_new,grid_data ,  'edgecolor','none')
                    %caxis([0,max(max(A0))]);
                    %surf(p.X_vol, p.Y_vol, A, 'edgecolor','none');
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
                    movegui(f2,[2350 700]);
                    clf(f2);
                    h6 = axes;
                    grid_data = griddata(p.X_vol,p.Y_vol, Rd,X_vol_new,Y_vol_new);
                    %surf(p.X_vol, p.Y_vol, Rd);
                    surf(X_vol_new,Y_vol_new,grid_data ,  'edgecolor','none')
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
                
                %%%% Nutrients in the Sediment %%%%
                if(Rs_vid)
                    f3 = figure(3);
                    movegui(f3,[1750 170]);
                    clf(f3);
                    %h6 = axes;
                    plot(sed_points,Rs);
                    grid off
                    ylabel("Algae concentration [mgP/m^2]");
                    xlabel("distance from lake center [m]");
                    title("Sediment nutrient concentration [mgP/m^2]");
                    %set(h6, 'Ydir', 'reverse')
                    %axis([0 p.W 0 p.Lmax 0 max(max(A0))]) % Fixing window size for video
                    %axis([0 p.W 0 p.Lmax -2 2]) % Fixing window size for video
                    drawnow
                    %pause(0.0001)
                    frame = getframe(gcf);
                    if(Rs_vid) writeVideo(Rs_video,frame); end
                end
                
                
            end
            if(A_vid)  close(A_video);  end
            if(Rd_vid) close(Rd_video); end
            if(Rs_vid) close(Rs_video); end
        end
    end
end