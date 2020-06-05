
function [] = plot_results_fn_V4(t,Y_t,p)
% plots simulation results in one figure and saves them as a jpg.

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

% creation of refined grid for interpolation.
[X_vol_new,Y_vol_new] = grid_interpolation_fn(p,100,100);

%% Calculation of the light intensity
A_temp = A';
D_temp = D';
Z_vol_temp = p.Z_vol';
I_alt = p.I0 .* exp(-p.I_matrix*(p.kA.*A_temp(:) + p.kD.*D_temp(:)) - p.kbg.*Z_vol_temp(:));
I = reshape(I_alt,[(p.Xn-1) (p.Zn-1)])';
p.I = I;

%% Figure index
figureIndex = 1;
close all

%% Plot of simulation results

if(true)
    font_size = 18;
    
    %% individual plots
    % plot of phytoplankton nutrient density
    if(false)
        % Plots algal density where the values have been interpolated to
        % produce a smooth figure.
        grid_data_A = griddata(p.X_vol,p.Z_vol, p.q.*A,X_vol_new,Y_vol_new, 'v4'); % interpolates to refined grid.
        figure(figureIndex)
       % movegui(figureIndex,[1720 700]);
        figureIndex = figureIndex +1;
        h1 = axes;
        
        surf(X_vol_new,Y_vol_new,grid_data_A ,  'edgecolor','none') % cell
        caxis([0 inf]);
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(h1, 'Ydir', 'reverse')
        title("Phytoplankton [mgP/m^3]");
        ylabel("Depth [m]");
        xlabel("distance from lake center [m]");
        zlabel("Algae concentration [mgP/m^3]");
        colorbar
        colormap(jet)
        set(gca,'fontSize',font_size);
        z_theory_algae = 1/p.kbg *log((p.I0*(p.Gmax/(p.lbg_A+p.Ad) -1))/p.H); % theoretical maximum depth for benthic growth if light limited,
        z_theory_algae = round( z_theory_algae,2);
        % at this depth the light limited growth is equal to respiration losses and death rate.
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[max(max(A))+1, max(max(A))+1], '--','color','black','LineWidth',2); % plotting the max depth line
            mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            %uicontrol('Style','text','String',mystr,'fontsize',22)
            h = annotation('textbox','String',mystr,'fontsize',font_size-2,'EdgeColor','none');
            h.Position=[.46 .36 .3581 .0879];
            hold off
        end
        
    end
    
    % plot of dissolved nutrient density
    if(false)
        grid_data_Rd = griddata(p.X_vol,p.Z_vol, Rd,X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        figure(figureIndex)
        movegui(figureIndex,[2300 700]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,grid_data_Rd ,  'edgecolor','none')
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        set(gca,'fontSize',font_size);
        caxis([0 inf]);
        colorbar
        colormap(jet)
        title("Dissolved nutrients [mgP/m^3]");
        ylabel("Depth [m]");
        xlabel("distance from lake center [m]");
        
    end
    
    % plot of detritus nutrient density
    if(false)
        grid_data_D = griddata(p.X_vol,p.Z_vol, D, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        figure(figureIndex)
        movegui(figureIndex,[2880 700]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,grid_data_D ,  'edgecolor','none')
        %surf(p.X,p.Z,D_new,'edgecolor','none'); % manually interpolated values to grid nodes.
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        caxis([0 inf]);
        colorbar
        colormap(jet)
        title("Detritus [mgP/m^3]");
        ylabel("Depth [m]");
        xlabel("distance from lake center [m]");
        set(gca,'fontSize',font_size);
    end
    
    % plot of sediment nutrient density
    if(false)
        % creating points for sediment plot
        sed_points = zeros(1,p.Xn-1);
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %sed_interp = interp1(sed_points,Rs, query_points,'makima');
        sed_interp = interp1(sed_points,Rs, query_points,'pchip');
        figure(figureIndex)
        movegui(figureIndex,[1720 170]);
        figureIndex = figureIndex +1;
        plot( query_points,sed_interp);
        title("Sediment Nutrients [mgP/m^2]");
        ylabel("concentration [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        set(gca,'fontSize',font_size);
    end
    
    % plot of Benthic algal density [mgC/m^2]
    if(false)
        % creating points for benthic algae plot
        benth_points = zeros(1,p.Xn-1);
        benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %benth_interp = interp1(sed_points,B, query_points,'makima');
        benth_interp = interp1(sed_points,B.*p.q_benth, query_points,'pchip'); % interpolated net growtht values, unit [mgP/(m^2 day)]
        figure(figureIndex)
        movegui(figureIndex,[2300 170]);
        figureIndex = figureIndex +1;
        x = [1:(p.W-1)/(p.Xn-1):p.W];
        plot(query_points, benth_interp);
        title("Benthic Algae [mgP/m^2]");
        ylabel("concentration [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        set(gca,'fontSize',font_size);
        hold on
        
        max_survival_depth = -1./p.kbg .* log(p.H_benth./( p.I0.*(p.Gmax_benth./p.lbg_benth -1) )); % maximum survival depth of benthic algae assuming no phytoplankton
        max_survival_depth_x = p.W*((max_survival_depth./p.Lmax -1)*(p.Lmax/(p.Lmin-p.Lmax)))^(1/p.alpha); % corressponding distance from the lake center
        
        
        if(  max_survival_depth_x  < p.Lmax && max_survival_depth_x > 0)
            xline(max_survival_depth_x, '--'); % marking out the maximum survival depth as a vertical line.
        end
        
        % plot of maximal benthic concentration at each point given light
        % limitation only.
        if(true)
            % Available light at the very bottom of the lake, at the center of the
            % bottom border of each cell. (previously the light at the center of
            % each bottom cell was used, but the light is attenuated a bit more
            % than that at the very bottom.)
            I_bottom = p.I(end,:).*exp(-( 0.5*(p.Z(end,1:end-1)+ p.Z(end,2:end)) - p.Z_vol(end,:)).*(p.kA.*A(end,:) + p.kD.*D(end,:) + p.kbg));
            
            syms b
            B_max = zeros(1, p.Xn-1); % maximum benthic concentration at each point if light limited
            for i =1:p.Xn-1
                eq =   p.Gmax_benth./p.kB.*log((p.H_benth + I_bottom(i))./(p.H_benth + I_bottom(i).*exp(-p.kB.*b))) == p.lbg_benth.*b ;
                B_max(i) =  vpasolve(eq,b,1e20);
            end
            max_benth_interp = interp1(sed_points,B_max.*p.q_benth, query_points,'pchip');
            
            % plot(query_points, max_benth_interp);
            hold off
        end
        
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
        ylabel("Depth [m]");
        xlabel("distance from lake center [m]");
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
    if(false)
        
        % calculation  of light/nutrient limitation at each point in lake.
        
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        
        interp_l_lim = griddata(p.X_vol,p.Z_vol, l_lim, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        interp_n_lim = griddata(p.X_vol,p.Z_vol, n_lim, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        
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
        
        z_theory_algae = 1/p.kbg *log((p.I0*(p.Gmax/(p.lbg_A+p.Ad) -1))/p.H); % theoretical maximum depth for benthic growth if light limited,
        % at this depth the light limited growth is equal to respiration losses and death rate.
        z_theory_algae = round( z_theory_algae,2);
        
        figure(figureIndex)
        movegui(figureIndex,[2880 170]);
        figureIndex = figureIndex +1;
        surf(X_vol_new,Y_vol_new,limitation ,  'edgecolor','none')
        grid off
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colormap(jet)
        title("light vs nutrient limitation");
        xlabel("red = nutrient limited, blue = light limited.");
        ylabel("Depth [m]");
        
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[2 2], '--','color','black','LineWidth',2);
            mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            h = annotation('textbox','String',mystr,'fontsize',font_size -2,'EdgeColor','none');
            h.Position=[.54 .37 .3581 .0879];
            hold off
            set(gca,'fontSize',font_size);
        end
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
    
    % plot of nutrient leakage over time
    if(false)
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
            leakage(i)      = (p.ntot_0 - ntot_end)/p.ntot_0;
        end
        figure(figureIndex)
        figureIndex = figureIndex +1;
        plot(t,leakage);
        %plot(leakage);
        title('(ntot(t)-ntot0)/ ntot0');
    end
    
    % plot of benthic growth
    if(false)
        
        figure(figureIndex)
        figureIndex = figureIndex +1;
        movegui(gca,[580 170]);
        clf(gca);
        
        % Bottom sediment interaction and benthic algal growth
        i = p.Zn-1; % eta
        j = 1:p.Xn-1; % xi
        
        % Available light at the very bottom of the lake, at the center of the
        % bottom border of each cell. (previously the light at the center of
        % each bottom cell was used, but the light is attenuated a bit more
        % than that at the very bottom.)
        B = B';
        I_bottom = p.I(end,:).*exp(-( 0.5*(p.Z(end,1:end-1)+ p.Z(end,2:end)) - p.Z_vol(end,:)).*(p.kA.*A(end,:) + p.kD.*D(end,:) + p.kbg));
        
        % Bethic algae net growth
        nutrient_limited_growth =  p.Gmax_benth.*B(j)'.*(Rd(i,j)./(Rd(i,j)+p.M_benth));
        light_limited_growth = p.Gmax_benth./p.kB.*log((p.H_benth + I_bottom)./(p.H_benth + I_bottom.*exp(-p.kB.*B(j)')));
        
        dBdt = zeros(1,p.Xn-1);
        dBdt(j) =  min(nutrient_limited_growth, light_limited_growth) - p.lbg_benth*B(j)';
        
        % creating points for benthic algae plot
        benth_points = zeros(1,p.Xn-1);
        benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %benth_interp = interp1(sed_points,B, query_points,'makima');
        net_growth_interp = interp1(sed_points,p.q_benth.*dBdt,query_points,'pchip'); % net benthic growth in [mgp/m^2 day]
        light_limited_interp =  interp1(sed_points,p.q_benth.*(light_limited_growth - p.lbg_benth*B(j)'),query_points,'pchip') ; % net benthic growth in [mgp/m^2 day] if light limited only
        nutrient_limited_interp =  interp1(sed_points,p.q_benth.*( nutrient_limited_growth - p.lbg_benth*B(j)'),query_points,'pchip'); % net benthic growth in [mgp/m^2 day]  if nutrient limited only
        plot(query_points,net_growth_interp,query_points,light_limited_interp,query_points,nutrient_limited_interp);
        grid off
        ylabel("growth [mgP/(m^2 day)]");
        xlabel("distance from lake center [m]");
        title("Benthic growth");
        ylim([0 inf]);
        legend("net growth", "light limited", "nutrient limited*");
        set(gca,'FontSize',font_size)
    end
    
    %% plots of all concentrations in one figure using griddata
    if(true)
        figure(figureIndex)
        figureIndex = figureIndex +1;
        set(gcf, 'Position',  [1730, 50, 1700, 1300]) % for widescreen
        % set(gcf, 'Position',  [0, 0, 1200, 1000]) % if on laptop
        tile_fig =  tiledlayout(3,3,'TileSpacing','Compact');
        tile_fig.Padding = 'compact';
        title_str = "";
        %title_str = "dx: " + num2str(p.dx(1,1)) + "  dz: " + num2str(p.dz(1,1)) + "  Xn: " + num2str(p.Xn) + "  Zn: " + num2str(p.Zn) + "  alpha: " + num2str(p.alpha) + "  kbg: " + num2str(p.kbg) + "";
        if(exist('p.stratified'))
            if(p.stratified)
                title_str = "Stratified, ";
            end
        end
        title_str = title_str + "Xn: " + num2str(p.Xn) + "  Zn: " + num2str(p.Zn) + "  alpha: " + num2str(p.alpha) + "  kbg: " + num2str(p.kbg) + "";
        title(tile_fig,title_str,'FontSize', font_size + 12); % shared title for all plots
        x_aspect = 1.2;
        
        
        %% phytoplankton
        
        % Plots algal density where the values have been interpolated to
        % produce a smooth figure.
        grid_data_A = griddata(p.X_vol,p.Z_vol, p.q.*A,X_vol_new,Y_vol_new, 'v4'); % interpolates to refined grid.
        nexttile
        surf(X_vol_new,Y_vol_new,grid_data_A ,  'edgecolor','none') % cell
        caxis([0 3.3]);
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        
        set(gca, 'Ydir', 'reverse')
        title("Phytoplankton [mgP/m^3]");
        ylabel("Depth [m]");
        %xlabel("distance from lake center [m]");
        zlabel("Algae concentration [mgP/m^3]");
        colorbar
        colormap(jet)
        pbaspect([1.2 1 1])
        set(gca,'fontSize',font_size);
        z_theory_algae = 1/p.kbg *log((p.I0*(p.Gmax/(p.lbg_A+p.Ad) -1))/p.H); % theoretical maximum depth for benthic growth if light limited,
        z_theory_algae = round( z_theory_algae,2);
        % at this depth the light limited growth is equal to respiration losses and death rate.
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[max(max(A))+1, max(max(A))+1], '--','color','black','LineWidth',2); % plotting the max depth line
            mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            h = annotation('textbox','String',mystr,'fontsize',font_size -2 ,'EdgeColor','none');
            h.Position=[.115 .635 .3581 .0879];
            hold off
        end
        
        %% Dissolved nutrients
        nexttile
        grid_data_Rd = griddata(p.X_vol,p.Z_vol, Rd,X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        surf(X_vol_new,Y_vol_new,grid_data_Rd ,  'edgecolor','none')
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        pbaspect([x_aspect 1 1])
        set(gca, 'Ydir', 'reverse')
        set(gca,'fontSize',font_size);
        caxis([0 inf]);
        colorbar
        colormap(jet)
        title("Dissolved nutrients [mgP/m^3]");
        %ylabel("Depth [m]");
        %xlabel("distance from lake center [m]");
        
        %% Detritus
        grid_data_D = griddata(p.X_vol,p.Z_vol, D, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        nexttile
        surf(X_vol_new,Y_vol_new,grid_data_D ,  'edgecolor','none')
        %surf(p.X,p.Z,D_new,'edgecolor','none'); % manually interpolated values to grid nodes.
        shading interp
        grid off
        az =0;
        el = 90;
        view(az,el);
        pbaspect([x_aspect 1 1])
        set(gca, 'Ydir', 'reverse')
        caxis([0 inf]);
        colorbar
        colormap(jet)
        title("Detritus [mgP/m^3]");
        % ylabel("Depth [m]");
        % xlabel("distance from lake center [m]");
        set(gca,'fontSize',font_size);
        
        %% Sediment nutrients
        sed_points = zeros(1,p.Xn-1); % creating points for sediment plot
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %sed_interp = interp1(sed_points,Rs, query_points,'makima');
        sed_interp = interp1(sed_points,Rs, query_points,'pchip');
        nexttile
        plot( query_points,sed_interp);
        pbaspect([x_aspect 1 1])
        title("Sediment Nutrients [mgP/m^2]");
        ylabel("concentration [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        set(gca,'fontSize',font_size);
        
        %% Benthic algae
        % creating points for benthic algae plot
        benth_points = zeros(1,p.Xn-1);
        benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %benth_interp = interp1(sed_points,B, query_points,'makima');
        benth_interp = interp1(sed_points,B.*p.q_benth, query_points,'pchip'); % interpolated net growtht values, unit [mgP/(m^2 day)]
        x = [1:(p.W-1)/(p.Xn-1):p.W];
        nexttile
        plot(query_points, benth_interp);
        pbaspect([x_aspect 1 1])
        title("Benthic Algae [mgP/m^2]");
        ylabel("concentration [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        set(gca,'fontSize',font_size);
        hold on
        
        max_survival_depth = -1./p.kbg .* log(p.H_benth./( p.I0.*(p.Gmax_benth./p.lbg_benth -1) )); % maximum survival depth of benthic algae assuming no phytoplankton
        max_survival_depth_x = p.W*((max_survival_depth./p.Lmax -1)*(p.Lmax/(p.Lmin-p.Lmax)))^(1/p.alpha); % corressponding distance from the lake center
        
        
        if(  max_survival_depth_x  < p.Lmax && max_survival_depth_x > 0)
            xline(max_survival_depth_x, '--'); % marking out the maximum survival depth as a vertical line.
        end
        
        % plot of maximal benthic concentration at each point given light
        % limitation only.
        if(true)
            % Available light at the very bottom of the lake, at the center of the
            % bottom border of each cell. (previously the light at the center of
            % each bottom cell was used, but the light is attenuated a bit more
            % than that at the very bottom.)
            I_bottom = p.I(end,:).*exp(-( 0.5*(p.Z(end,1:end-1)+ p.Z(end,2:end)) - p.Z_vol(end,:)).*(p.kA.*A(end,:) + p.kD.*D(end,:) + p.kbg));
            
            syms b
            B_max = zeros(1, p.Xn-1); % maximum benthic concentration at each point if light limited
            for i =1:p.Xn-1
                eq =   p.Gmax_benth./p.kB.*log((p.H_benth + I_bottom(i))./(p.H_benth + I_bottom(i).*exp(-p.kB.*b))) == p.lbg_benth.*b ;
                B_max(i) =  vpasolve(eq,b,1e20);
            end
            max_benth_interp = interp1(sed_points,B_max.*p.q_benth, query_points,'pchip');
            
            % plot(query_points, max_benth_interp);
            hold off
        end
        
        %% light/nutrient limitation
        % calculation  of light/nutrient limitation at each point in lake.
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        interp_l_lim = griddata(p.X_vol,p.Z_vol, l_lim, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        interp_n_lim = griddata(p.X_vol,p.Z_vol, n_lim, X_vol_new,Y_vol_new,'v4'); % interpolates to refined grid.
        
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
        
        z_theory_algae = 1/p.kbg *log((p.I0*(p.Gmax/(p.lbg_A+p.Ad) -1))/p.H); % theoretical maximum depth for benthic growth if light limited,
        % at this depth the light limited growth is equal to respiration losses and death rate.
        z_theory_algae = round( z_theory_algae,2);
        
        nexttile
        surf(X_vol_new,Y_vol_new,limitation ,  'edgecolor','none')
        grid off
        pbaspect([x_aspect 1 1])
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colormap(jet)
        title("light vs nutrient limitation");
        xlabel("distance from lake center [m]");
        ylabel("Depth [m]");
        set(gca,'fontSize',font_size);
        map = [0 0 0.75; 0.75 0 0];
        colormap(gca,map);
        
        dim = [.82 .15 .3 .3];
        str = '   red - nutrient limited \n blue - light limited';
        annotation('textbox',dim,'String',sprintf(str),'FitBoxToText','on','EdgeColor','none','FontSize',font_size-2);
        
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[2 2], '--','color','black','LineWidth',2);
            %mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            %h = annotation('textbox','String',mystr,'fontsize',15,'EdgeColor','none');
            %h.Position=[.54 .37 .3581 .0879];
            hold off
            set(gca,'fontSize',font_size);
        end
        
        %% nutrients proportions and losses
        t =  nexttile;
        delete(t);
        dim = [.11 .1 .2 .2];
        perc_algae = sum(sum(n_algae_end))/ntot_end;
        perc_dissolved = sum(sum(n_dissolved_end))/ntot_end;
        perc_detritus = sum(sum(n_detritus_end))/ntot_end;
        perc_sediment = sum(n_sediment_end)/ntot_end;
        perc_benthic = sum(n_benthic_end)/ntot_end;
        [x,isterm,dir] = eventfun_V4(t(end),Y_t(end,:)',p);
        str = '';
        str  = append( '______________________________________________', '\n');
        str  = append(str, 'Norm of the time derivatives at steady state: ', num2str(x), ' \n ');
        str  = append(str , 'portion of leaked nutrients: ', num2str(abs((p.ntot_0-ntot_end) /p.ntot_0 )), ' \n ' );
        str  = append(str ,'leaked nutrients [mgP]:      ', num2str(abs(p.ntot_0-ntot_end)), ' \n ');
        str  = append(str , '______________________________________________', '\n');
        str  = append(str ,'portion of nutrients in phytoplankton:        ', num2str(perc_algae), ' \n ');
        str  = append(str ,'portion of nutrients in dissolved nutrients: ', num2str(perc_dissolved), ' \n ');
        str  = append(str ,'portion of nutrients in detritus:                   ', num2str(perc_detritus), ' \n ');
        str  = append(str ,'portion of nutrients in sediment:                ', num2str(perc_sediment), ' \n ');
        str  = append(str ,'portion of nutrients in benthic algae:         ', num2str(perc_benthic));
        textbox = annotation('textbox',dim,'String',sprintf(str),'FitBoxToText','on','EdgeColor','none');
        textbox.FontSize = font_size-4;
    end
    
    %% save figure
    alpha_str = "VALUE";
    if(p.alpha == 1.5) alpha_str = "1_5"; end
    if(p.alpha == 1) alpha_str = "1"; end
    
    kbg_str = "VALUE";
    if(p.kbg == 0) kbg_str = "0"; end
    if(p.kbg == 0.02)kbg_str = "0_02"; end
    if(p.kbg == 0.2)kbg_str = "0_2"; end
    if(p.kbg == 0.4)kbg_str = "0_4"; end
    if(p.kbg == 0.8)kbg_str = "0_8"; end
    if(p.kbg == 2.0)kbg_str = "2"; end
    
    if(p.resus == 0) resus_str = "0_0"; end
    if(p.resus == 0.1) resus_str = "0_1"; end
    if(p.resus == 0.5) resus_str = "0_5"; end
    
    if(p.r == 0) remin_str = "0_0"; end
    if(p.r == 0.05) remin_str = "0_05"; end
    if(p.r == 0.1) remin_str = "0_1"; end
    
    
    %fig_name = "2D_benthic_results_V4_dx_" + num2str(p.dx(1,1)) + "_dz_" + num2str(p.dz(1,1)) + "_Xn_" + num2str(p.Xn) + "_Zn_" + num2str(p.Zn) + "_alpha_" + alpha_str + "_kbg_"+ kbg_str + "";
    %fig_name = "2D_benthic_results_V4_dx_" + num2str(p.dx(1,1)) + "_dz_" + num2str(p.dz(1,1)) + "_Xn_" + num2str(p.Xn) + "_Zn_" + num2str(p.Zn) + "_alpha_" + alpha_str + "_kbg_"+ kbg_str + "_resuspRate_" + resus_str + "";
    %fig_name = "2D_benthic_results_V4_dx_" +num2str(p.dx(1,1)) + "_dz_" + num2str(p.dz(1,1)) + "_Xn_" + num2str(p.Xn) + "_Zn_" + num2str(p.Zn) + "_alpha_" + alpha_str + "_kbg_"+ kbg_str +  "_reminrate_"+ remin_str ;
     fig_name = "2D_benthic_results_V4_dx_" + "Xn_" + num2str(p.Xn) + "_Zn_" + num2str(p.Zn) + "_alpha_" + alpha_str + "_kbg_"+ kbg_str +"_therm_depth_10_therm_thickness_3";
     
    saveas(gcf,fig_name,'jpg');
end
end