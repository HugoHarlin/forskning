%  clc
%  clear
%  close all

% Plots the nutrient of each vertical water column as a heat map.
% Algae, dissolved nutrients, and the sediment is plotted separately as
% well as the sum of the three.

%% Loads workspace and generates nutrient plots

 %dx=10000;
 %dy=0.5;
% load("dx_"+ num2str(dx)+ "_dy_" + num2str(dy) + "_res_50_50");

N_A = sum(n_algae_end);
N_R = sum(n_dissolved_end);
N_Rs = n_sediment_end;
sh = 0;
fontsize = 18;
p.dx
p.dy


    
figure(1)
x0=350;
y0=120;
width=550;
height=450;
set(gcf,'position',[x0,y0,width,height])


%% heatmap style plots
if(false)
    %% Algae
    subaxis(4, 1, 1, 'sh', sh, 'sv', 0.00, ...
        'PL', 0.10,  'PR', 0.20, 'PT', 0.01, 'PB', 0.01, ...
        'MT', 0.0,'MB', 0.080, 'ML', 0, 'MR', 0.0);
    imagesc(N_A);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca, 'Xdir', 'reverse')
    set(gca,'FontSize', fontsize);
    colormap(jet);
    cb1 = colorbar;
    colormap(jet);
%     lower = cb1.Limits(1) + 0.195*(mean(cb1.Limits));
%     upper = cb1.Limits(2) - 0.195*(mean(cb1.Limits));
%     set(cb1, 'YTick', round(linspace(lower, upper, 3)));
    %set(cb1, 'YTick', linspace(lower, upper, 3));
    % set(cb1, 'YTick', [0.2 0.8 1.5]);
    
    %% Dissolved nutrients
    subaxis(4, 1, 2, 'sh', sh, 'sv', 0.00, ...
        'PL', 0.10,  'PR', 0.20, 'PT', 0.01, 'PB', 0.01, ...
        'MT', 0.0,'MB', 0.080, 'ML', 0, 'MR', 0.0);
    colormap(jet);
    imagesc(N_R);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca, 'Xdir', 'reverse')
    set(gca,'FontSize', fontsize);
    
    cb2 = colorbar;
    colormap(jet);
%     lower = cb2.Limits(1) + 0.15*(mean(cb2.Limits));
%     upper = cb2.Limits(2) - 0.15*(mean(cb2.Limits));
%     set(cb2, 'YTick', round(linspace(lower, upper, 3)));
    
    %% Sediment
    subaxis(4, 1, 3, 'sh', sh, 'sv', 0.00, ...
        'PL', 0.10,  'PR', 0.20, 'PT', 0.01, 'PB', 0.03, ...
        'MT', 0.0,'MB', 0.080, 'ML', 0, 'MR', 0.0);
    imagesc(N_Rs);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    set(gca, 'Xdir', 'reverse')
    set(gca,'FontSize', fontsize);
    colormap(jet);
    cb3 = colorbar;
    colormap(jet);
%     lower = cb3.Limits(1) + 0.0085*(mean(cb3.Limits));
%     upper = cb3.Limits(2) - 0.0085*(mean(cb3.Limits));
%     set(cb3, 'YTick', round(linspace(lower, upper, 3)));
%     
    %% total
    subaxis(4, 1, 4, 'sh', sh, 'sv', 0.00, ...
        'PL', 0.10,  'PR', 0.20, 'PT', 0.01, 'PB', 0.02, ...
        'MT', 0.0,'MB', 0.080, 'ML', 0, 'MR', 0);
    imagesc( N_A + N_R + N_Rs);
    set(gca,'YTick',[]);
    set(gca, 'Xdir', 'reverse')
    if(dx~=1000)
        set(gca,'XTick',[]);
    end
    
    set(gca,'FontSize', fontsize);
    colormap(jet);
    cb4 = colorbar;
    colormap(jet);
%     lower = cb4.Limits(1) + 0.0085*(mean(cb4.Limits));
%     upper = cb4.Limits(2) - 0.0085*(mean(cb4.Limits));
%     set(cb4, 'YTick', round(linspace(lower, upper, 3)));
%     
end

%% Log scale plots
if(true)
x = linspace(1,50,50);
N_A = fliplr(N_A);
N_R = fliplr(N_R);
N_Rs = fliplr(N_Rs);

plot(x,N_A,x,N_R,x,N_Rs, 'linewidth', 2);
set(gca, 'YScale', 'log');
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 1);
xlim([1 50]) 
ylim([1,1E4])
grid on
lgd = legend('Algae','Dissolved Nutrients','Sediment');
lgd.FontSize = 13;


if(dx ~= 1000)
    set(gca,'XTickLabel',[]);
end

disp("algal nutrients / total nutrients: " + num2str(sum(N_A) / ntot_0));
disp("dissolved nutrients / total nutrients: " + num2str(sum(N_R) / ntot_0));
disp("sediment nutrients / total nutrients: " + num2str(sum(N_Rs) / ntot_0));
end