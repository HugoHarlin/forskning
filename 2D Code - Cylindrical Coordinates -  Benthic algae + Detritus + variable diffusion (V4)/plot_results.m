
%% Plots the nutrient contents of the state variables

%% Loads workspace and generates nutrient plots
% clc
clear
close all

%% loop over the diffusion coefficients
res = 26;
dx=10;
dz=100;

load("2D_benthic_results_V4_dx_" + num2str(dx) + "_dz_" +num2str(dz)+ "_Xn_" +num2str(res)+"_Zn_"+ num2str(res));

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

%% heatmap of nutrinet percentages
% matrixes for storing the nutrient content percentage

perc_algae_m = zeros(4);
perc_dissolved_m = zeros(4);
perc_detritus_m = zeros(4);
perc_sediment_m = zeros(4);
perc_benthic_m = zeros(4);
res = 26;

for x = 1 :4
    for z = 1:4
        
        dx = 10^(x-1);
        dz = 10^(z-1);
        
        load("2D_benthic_results_V4_dx_" + num2str(dx) + "_dz_" +num2str(dz)+ "_Xn_" +num2str(res)+"_Zn_"+ num2str(res));
        
        % extracting equilibrium concentrations
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
        
        % calculating nutrient percentages
        perc_algae = sum(sum(n_algae_end))/ntot_end;
        perc_dissolved = sum(sum(n_dissolved_end))/ntot_end;
        perc_detritus = sum(sum(n_detritus_end))/ntot_end;
        perc_sediment = sum(n_sediment_end)/ntot_end;
        perc_benthic = sum(n_benthic_end)/ntot_end;
        
        % storing results
        perc_algae_m(z,x) = perc_algae;
        perc_dissolved_m(z,x) = perc_dissolved;
        perc_detritus_m(z,x) = perc_detritus;
        perc_sediment_m(z,x) = perc_sediment;
        perc_benthic_m(z,x) = perc_benthic;
        
    end
end


%% heatmaps of the nutrient percentages
figure(1)
subplot(2,3,1) % algae
imagesc([1,2,3,4],[1,2,3,4],perc_algae_m);
xticks([1,2,3,4]);
yticks([1,2,3,4]);
xticklabels({'1','10','100','1000'});
yticklabels({'1','10','100','1000'});
colorbar
xlabel('dx [m^2/day]');
ylabel('dz [m^2/day]');
title('Phytoplankton');
%caxis([0 0.6]);

subplot(2,3,2) % dissolved nutrients
imagesc([1,2,3,4],[1,2,3,4],perc_dissolved_m);
xticks([1,2,3,4]);
yticks([1,2,3,4]);
xticklabels({'1','10','100','1000'});
yticklabels({'1','10','100','1000'});
colorbar
xlabel('dx [m^2/day]');
ylabel('dz [m^2/day]');
title('Dissolved nutrients');
%caxis([0 0.6]);

subplot(2,3,3) % detritus
imagesc([1,2,3,4],[1,2,3,4],perc_detritus_m);
xticks([1,2,3,4]);
yticks([1,2,3,4]);
xticklabels({'1','10','100','1000'});
yticklabels({'1','10','100','1000'});
colorbar
xlabel('dx [m^2/day]');
ylabel('dz [m^2/day]');
title('Detritus');
%caxis([0 0.6]);

subplot(2,3,4) % sediment
imagesc([1,2,3,4],[1,2,3,4],perc_sediment_m);
xticks([1,2,3,4]);
yticks([1,2,3,4]);
xticklabels({'1','10','100','1000'});
yticklabels({'1','10','100','1000'});
colorbar
xlabel('dx [m^2/day]');
ylabel('dz [m^2/day]');
title('Sediment');
%caxis([0 0.6]);

subplot(2,3,5) % benthic algae
imagesc([1,2,3,4],[1,2,3,4],perc_benthic_m);
xticks([1,2,3,4]);
yticks([1,2,3,4]);
xticklabels({'1','10','100','1000'});
yticklabels({'1','10','100','1000'});
colorbar
xlabel('dx [m^2/day]');
ylabel('dz [m^2/day]');
title('Benthic Algae');
%caxis([0 0.6]);

