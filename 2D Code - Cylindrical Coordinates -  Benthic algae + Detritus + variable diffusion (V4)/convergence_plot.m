
%%
clc
clear
close all

%res_vec = [5:26,36,41]; %dx 10 dz 10
%res_vec = [5:22,33]; % dx 100 dz 100
res_vec = [16,19,23,26,31]; %,35]; % straatified kbg 0.8 therm depth 10 m

num_res = length(res_vec);
max_res = max(res_vec);
detritus_tot = zeros(1,num_res); % totalt nutrient content in detritus for the simulated resolutions
algae_tot = zeros(1,num_res); % percentage of total nutrients in algae
dissolved_tot = zeros(1,num_res); % percentage of total nutrients in dissolved state
sed_tot = zeros(1,num_res); % percentage of total nutrients in sediment
benth_tot = zeros(1,num_res); % percentage of total nutrients in benthic algae

detr_matrix = zeros(num_res,max_res,max_res); % detritus nutrient concentration at each point for all resolutions 
nutr_matrix = zeros(num_res,max_res,max_res); % dissolved nutrient concentration at each point for all resolutions 
algae_matrix = zeros(num_res,max_res,max_res); % phytoplankton nutrient concentration at each point for all resolutions 
benth_vec = zeros(num_res,max_res);  % benthic algae nutrient concentration at each point for all resolutions 
sed_vec = zeros(num_res,max_res);  % sediment nutrient concentration at each point for all resolutions 

for i = 1:num_res % interating over the resolutions
    
    % declaring struct
    eval("p_" + num2str(res_vec(i)) + "= struct;");
        
    % loading data into struct
    % eval("p_" + num2str(res_vec(i)) + "=
    % load(""2D_convergence_res_V4_Xn_" + num2str(res_vec(i)) + "_Zn_" + num2str(res_vec(i))+ """);"); 
    eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V4_Xn_" + num2str(res_vec(i)) + "_Zn_" + num2str(res_vec(i))+ "_alpha_1_5_kbg_0_8_Stratified_therm_depth_10"");"); % stratified simulations
    
    % Extracting  algae, dissolved nutrients, and detritus concentrations
    eval(" A = p_" + num2str(res_vec(i)) + ".Y_t(end,(1:(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1)));");
    eval("Rd = p_" + num2str(res_vec(i)) + ".Y_t(end,(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1)+1 : 2*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1));");
    eval(" D = p_" + num2str(res_vec(i)) + ".Y_t(end,2*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1)+1 : 3*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1));");
    eval("Rs = p_" + num2str(res_vec(i)) + ".Y_t(end,3*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1)+1 : 3*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1) +p_" + num2str(res_vec(i)) + ".p.Xn-1 );");
    eval(" B = p_" + num2str(res_vec(i)) + ".Y_t(end,3*(p_" + num2str(res_vec(i)) + ".p.Xn-1)*(p_" + num2str(res_vec(i)) + ".p.Zn-1)+ (p_" + num2str(res_vec(i)) + ".p.Xn-1) +1 : end);");
    
    
    % reformatting into a square matrix
    eval("A = reshape(A, [p_" + num2str(res_vec(i)) + ".p.Xn-1, p_" + num2str(res_vec(i)) + ".p.Zn-1]);");
    eval("D = reshape(D, [p_" + num2str(res_vec(i)) + ".p.Xn-1, p_" + num2str(res_vec(i)) + ".p.Zn-1]);");
    eval("Rd = reshape(Rd, [p_" + num2str(res_vec(i)) + ".p.Xn-1, p_" + num2str(res_vec(i)) + ".p.Zn-1]);");
    A = A';
    D = D';
    Rd = Rd';
    
    % calculating nutrient content for the state vectors in each grid
    % element.
    eval("algae = sum(sum( p_"+num2str(res_vec(i))+".p.volumes_cyl.*A.*p_"+ num2str(res_vec(i))+".p.q));");
    eval("detr = sum(sum( p_"+num2str(res_vec(i))+".p.volumes_cyl.*D));");
    eval("nutr = sum(sum( p_"+num2str(res_vec(i))+".p.volumes_cyl.*Rd));");
    eval("benth = sum(sum( p_"+num2str(res_vec(i))+".p.Area_bottom_cyl.*p_"+ num2str(res_vec(i))+".p.q.*B));");
    eval("sed = sum(sum( p_"+num2str(res_vec(i))+".p.Area_bottom_cyl.*Rs));");
    
    % storing total nutrient concentration for each state variable
    detritus_tot(i) = sum(sum(detr)); 
    dissolved_tot(i) = sum(sum(nutr));
    algae_tot(i) = sum(sum(algae));
    sed_tot(i) = sum(sed);
    benth_tot(i) = sum(benth);
    
    % Interpolating to a refined grid (same grid as the highest resolution)
    query_points = eval("linspace(0,p_" + num2str(res_vec(i)) + ".p.W,max_res);");
    sed_points = eval("zeros(1, p_"+num2str(res_vec(i))+".p.Xn-1);");
    sed_points(:) = eval("( p_"+num2str(res_vec(i))+".p.X(1,2:end) +  p_"+num2str(res_vec(i))+".p.X(1,1:end-1))/2;");
    [X_vol_new,Y_vol_new] = eval("grid_interpolation_fn(p_"+num2str(res_vec(i))+".p,max_res,max_res);"); % refined grid for interpolation.
    eval("algae_matrix("+num2str(i)+",:,:) = griddata(p_"+num2str(res_vec(i))+".p.X_vol, p_"+num2str(res_vec(i))+".p.Z_vol, A.*p_"+ num2str(res_vec(i))+".p.q,X_vol_new,Y_vol_new, 'v4');");
    eval("detr_matrix("+num2str(i)+",:,:) = griddata(p_"+num2str(res_vec(i))+".p.X_vol, p_"+num2str(res_vec(i))+".p.Z_vol, D,X_vol_new,Y_vol_new, 'v4');");
    eval("nutr_matrix("+num2str(i)+",:,:) = griddata(p_"+num2str(res_vec(i))+".p.X_vol, p_"+num2str(res_vec(i))+".p.Z_vol, Rd,X_vol_new,Y_vol_new, 'v4');");
    eval("sed_vec(i,:) = interp1(sed_points,Rs, query_points,'pchip');");
    eval("benth_vec(i,:) = interp1(sed_points,B, query_points,'pchip');");
end

%% %%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%
close all

% sediment nutrients
figure(1)
for i =5:max_res % i loops over each point in the sediment layer, i=1 corresponds to the leftmost point at the center of the lake.
    %plot(res_vec, abs((sed_vec(:,i)-sed_vec(end,i))./sed_vec(end,i)));
    plot(res_vec,sed_vec(:,i)/(sed_vec(end,i)));
    hold on
end
title("Local comparison of the sediment nutrient content (normalized).");
hold off

% benthic algae
figure(2)
for i =5:max_res % i loops over each point in the benthic layer,  i=1 corresponds to the leftmost point at the center of the lake.
    plot(res_vec, benth_vec(:,i)./benth_vec(end,i));
    hold on
end
title("Local comparison of the benthic nutrient content (normalized).");
hold off


% plots of local convergence detritus
figure(3)
for i =1:max_res
    for j =1:max_res
        plot(res_vec, detr_matrix(:,i,j)./detr_matrix(end,i,j));
        hold on
    end
end
title("Local comparison of the detritus nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of local convergence dissolved nutrients
figure(4)
for i =1:max_res
    for j =1:max_res
        plot(res_vec, nutr_matrix(:,i,j)./nutr_matrix(end,i,j));
        hold on
    end
end
title("Local comparison of the dissolved nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of local convergence dissolved nutrients
figure(5)
for i =1:max_res
    for j =1:max_res
        plot(res_vec, algae_matrix(:,i,j)./algae_matrix(end,i,j));
        hold on
    end
end
title("Local comparison of the plankton nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of the total nutrient content percentage of each state variable
figure(6)
%plot(res_vec, detritus_tot./detritus_tot(end), res_vec, dissolved_tot./dissolved_tot(end), res_vec, algae_tot./algae_tot(end));
plot(res_vec, detritus_tot./detritus_tot(end), res_vec, dissolved_tot./dissolved_tot(end));
hold on
plot(res_vec, benth_tot./benth_tot(end), res_vec, sed_tot./sed_tot(end));
%legend("detritus","dissolved nutrients","algae","benthic algae","sediment");
title("Normalized total nutrient of each state variable");
hold off



% function for interpolating resluts to a finer grid.
function [X_vol_new,Y_vol_new] = grid_interpolation_fn(p,res_x,res_z)
%% Mesh generator
%Generates finer mesh for interpolation of results, makes for nicer plots.
%
% Hugo Harlin 2019

% Quantities relating to system size
p.Xn_new = res_x+1; % Number of grid-points (width)
p.Yn_new = res_z+1; % Number of grid-points (depth)

% Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
% at the center of the lake ( slope = alpha* (Lmin - Lmax)/W  ).
% (0,0) is placed at the center of the lake at the surface, y-dim is facing
% downward and x-dim is facing towards the lake edge

p.X_new = zeros(p.Yn_new, p.Xn_new); % Mesh-spacing in x-dimension
p.Y_new = zeros(p.Yn_new, p.Xn_new); % Mesh-spacing in y-dimension

for i=1:1:p.Yn_new
    p.X_new(i,:) = [0:p.W/(p.Xn_new-1):p.W]; % even spacing of the grid in x-dimension
end

% the grid is compressed in y-dimension, with depth Lmin at the
% shore and Lmax at the center of the lake.
for i=1:p.Yn_new
    for j = 1:p.Xn_new
        p.Y_new(i,j) = (p.Lmax/(p.Yn_new-1))*(i-1)*(1 + (p.Lmin/p.Lmax -1)* (p.X_new(i,j)/p.W).^(p.alpha));
    end
end

% coordinates of the center of each mesh quadrilateral
% returns the coordinates of the center of each volume element in the mesh
X_vol_new = zeros(p.Yn_new-1, p.Xn_new-1);
Y_vol_new = zeros(p.Yn_new-1, p.Xn_new-1);


for j=1:p.Xn_new-1
    for i =1:p.Yn_new-1
        X_vol_new(i,j) = (p.X_new(i,j) + p.X_new(i+1,j) + p.X_new(i,j+1) + p.X_new(i+1,j+1) ) /4;
        Y_vol_new(i,j) = ( p.Y_new(i,j) +  p.Y_new(i+1,j) +  p.Y_new(i,j+1) +  p.Y_new(i+1,j+1) ) /4;
    end
end

end






