% clear all
% close all
% clc
% 
% % subplot of one set of figures
% 
dx=100000;
dy=10;
load("dx_"+ num2str(dx)+ "_dy_" + num2str(dy) + "_res_50_50");

    figure(1)
    x0=500;
    y0=100;
    width=550;
    height=850;
    set(gcf,'position',[x0,y0,width,height])
fontsize = 28;

xlabel = true;
ylabel = true;

%% Generating finer mesh for interpolation


% Quantities relating to system size
p.Xn = 201; % Number of grid-points (width)
p.Yn = 201; % Number of grid-points (depth)
p.Lmin = 0.1; % Minimum lake depth (depth at land-water interface) [m]
p.Lmax = 50; % Maximum lake depth [m]
p.W = 	50; % Lake width [m]
p.alpha = 1; % Exponent governing the slope of the lake bottom

% Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
% at the center of the lake ( slope = alpha* (Lmin - Lmax)/W  ).
% (0,0) is placed at the center of the lake at the surface, y-dim is facing
% downward and x-dim is facing towards the lake edge

p.X_new = zeros(p.Yn, p.Xn); % Mesh-spacing in x-dimension
p.Y_new = zeros(p.Yn, p.Xn); % Mesh-spacing in y-dimension

for i=1:1:p.Yn
    p.X_new(i,:) = [0:p.W/(p.Xn-1):p.W]; % even spacing of the grid in x-dimension
end

% the grid is compressed in y-dimension, with depth Lmin at the
% shore and Lmax at the center of the lake.
for i=1:p.Yn
    for j = 1:p.Xn
        p.Y_new(i,j) = (p.Lmax/(p.Yn-1))*(i-1)*(1 + (p.Lmin/p.Lmax -1)* (p.X_new(i,j)/p.W).^(p.alpha));
    end
end

% coordinates of the center of each mesh quadrilateral
% returns the coordinates of the center of each volume element in the mesh
X_vol_new = zeros(p.Yn-1, p.Xn-1);
Y_vol_new = zeros(p.Yn-1, p.Xn-1);


for j=1:p.Xn-1
    for i =1:p.Yn-1
        X_vol_new(i,j) = (p.X_new(i,j) + p.X_new(i+1,j) + p.X_new(i,j+1) + p.X_new(i+1,j+1) ) /4;
        Y_vol_new(i,j) = ( p.Y_new(i,j) +  p.Y_new(i+1,j) +  p.Y_new(i,j+1) +  p.Y_new(i+1,j+1) ) /4;
    end
end
p.X_vol_new = X_vol_new;
p.Y_vol_new = Y_vol_new;

%% Nutrients in sediment
sh = 0.05;

subaxis(7, 1, 1, 'sh', sh, 'sv', 0.00, ...
    'PL', 0.10,  'PR', 0.20, 'PT', 0.01, 'PB', 0.01, ...
    'MT', 0.0,'MB', 0.020, 'ML', 0, 'MR', 0.0);

imagesc(Rs);
 %imagesc(interp1(p.X_vol(1,:),Rs,p.X_vol_new(1,:)));

set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca, 'Xdir', 'reverse')
set(gca,'FontSize', fontsize);
%set(gca,'FontWeight', 'Bold');
colormap(jet);
cb = colorbar;
% lower = cb.Limits(1) + 0.00085*(mean(cb.Limits));
% upper = cb.Limits(2) - 0.00085*(mean(cb.Limits));
% set(cb, 'YTick', round(linspace(lower, upper, 3)));

%% Algae nutrient content

subaxis(7, 1, 2:4, 'sh', 0.000, 'sv', 0.000, ...
    'PL', 0.10,'PR', 0.2, 'PT', 0.000, 'PB', 0.02,...
    'MB', 0.030, 'MT', 0.000, 'ML', 0, 'MR', 0.00);

 surf(p.X_vol_new,p.Y_vol_new ,griddata(X_vol,Y_vol, p.q.*A,p.X_vol_new,p.Y_vol_new),  'edgecolor','none')
 %surf(p.X_vol,p.Y_vol,p.q.*A,  'edgecolor','none', 'FaceAlpha', 1)
view(az,el);
set(gca, 'YDir','reverse', 'XDir', 'reverse')
axis square
%set(gca, 'YDir','reverse')

% interval for common color bar
bottom = min(min(min(A)),min(min(Rd)));
top  = max(max(max(A)),max(max(Rd)));

% Common color bar
colormap(jet);
%         caxis manual
% caxis([bottom top]);


xticks(10:10:50);
yticks(10:10:50);
xlim([0 50])
ylim([0 50])

%%create a second axes.
ax1 = gca; % the first axes
ax1.FontSize = fontsize;
%ax1.FontWeight = 'bold';
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
ax2.FontSize = fontsize;
linkaxes([ax1 ax2],'xy')
grid on
axis square
temp = colorbar;
set(temp,'YTick',[])
set(gcf,'CurrentAxes',ax1); % return to ax1
cbar = colorbar;
% cbar.FontWeight = 'bold';
set(gca,'XTickLabel',[]);

if(dy~=1)
set(gca,'YTickLabel',[]);
end

%% Dissolved nutrients

subaxis(7, 1, 5:7, 'sh', 0.000, 'sv', 0.00, ...
    'PL', 0.100,'PR', 0.20, 'PT', 0.00, 'PB', 0.02,...
    'MT', 0.0,'MB', 0.05, 'ML', 0, 'MR', 0.00);

 surf(p.X_vol_new,p.Y_vol_new ,griddata(X_vol,Y_vol, Rd,p.X_vol_new,p.Y_vol_new),  'edgecolor','none')
alpha(1)
view(az,el);
%set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse', 'XDir', 'reverse')
view(2);
axis square
colormap(jet);
% set(gca, 'xtick' ,fliplr(get(gca, 'xtick')));

xticks(10:10:50);
yticks(10:10:50);
xlim([0 50])
ylim([0 50])
ax1 = gca; % the first axes
ax1.FontSize = fontsize;

%%create a second axes.
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

axis square
linkaxes([ax1 ax2],'xy')
grid on
temp = colorbar;
set(temp,'YTick',[])
set(gcf,'CurrentAxes',ax1); % return to ax1
cbar = colorbar;
cbar.FontSize = fontsize;

%Removed axis unless on the bottom/left edge
if(dx~=1000)
set(gca,'XTickLabel',[]);
end

if(dy~=1)
set(gca,'YTickLabel',[]);
end
