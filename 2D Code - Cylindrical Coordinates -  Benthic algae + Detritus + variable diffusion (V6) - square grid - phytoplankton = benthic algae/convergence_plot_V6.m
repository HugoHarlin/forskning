
%%
clc
clear
close all

% flags for choosing which parameter setup to plot
strafified_fixed_resus = 0;
strafified_fixed_resus_500m_radius = 1;
stratified_variable_resus_type_3_width_500 = 0;
stratified_variable_resus_type_1 = 0;
stratified_variable_resus_type_2 = 0;

folder = 'C:\Users\huha0015\Documents\PhD\Results V6\conv sim';
addpath(genpath(folder));
addpath('C:\Users\huha0015\Documents\PhD\Results V6\conv sim\conv sim variable resus\linear response function (type 2)');
addpath('C:\Users\huha0015\Documents\PhD\Results V6\conv sim\conv sim variable resus\quadratic response function (type 1)');
addpath('C:\Users\huha0015\Documents\PhD\Results V6\conv sim\conv sim fixed resus');
res_vec  = (10:4:58); % stratified, kbg = 0.2
%res_vec  = (12:2:34); % stratified, kbg = 0.2

num_res = length(res_vec);
max_res = max(res_vec);
detritus_tot = zeros(1,num_res); % totalt nutrient content in detritus for the simulated resolutions
algae_tot = zeros(1,num_res); % percentage of total nutrients in algae
dissolved_tot = zeros(1,num_res); % percentage of total nutrients in dissolved state
sed_tot = zeros(1,num_res); % percentage of total nutrients in sediment
benth_tot = zeros(1,num_res); % percentage of total nutrients in benthic algae

detr_matrix = zeros(num_res,sum(1:max_res-1)); % detritus nutrient concentration at each point for all resolutions
nutr_matrix = zeros(num_res,sum(1:max_res-1)); % dissolved nutrient concentration at each point for all resolutions
algae_matrix = zeros(num_res,sum(1:max_res-1)); % phytoplankton nutrient concentration at each point for all resolutions
benth_vec = zeros(num_res,max_res-1);  % benthic algae nutrient concentration at each point for all resolutions
sed_vec = zeros(num_res,max_res-1);  % sediment nutrient concentration at each point for all resolutions

ntot_vec = zeros(1,length(res_vec));

% coordinates of the grid centers for the highest resolution
% flattened coordinates of grid element centers
x_vol_vec_max = zeros(1,sum(1:max_res-1));
z_vol_vec_max = zeros(1,sum(1:max_res-1));

% load the highest resolution simulation into struct temp

if(strafified_fixed_resus)
    eval("temp = load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(max_res) +"_Stratified_therm_depth_2_thickness_3_profile_100_1_10_CONV_SIM"");"); % stratified simulations
end

if(stratified_variable_resus_type_1)
    eval("temp = load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(max_res) +"_Stratified_therm_depth_5_thickness_3_profile_100_1_10_variable_resus_resp_fn_1_CONVSIM"");"); % stratified simulations, variable resus
end

if(stratified_variable_resus_type_2)
    eval("temp = load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(max_res) +"_Stratified_therm_depth_5_thickness_3_profile_100_1_10_variable_resus_resp_fn_2_CONVSIM"");"); % stratified simulations, variable resus type 2
end

if(strafified_fixed_resus_500m_radius)
    eval("temp = load(""2D_benthic_results_V6_alpha_1_kbg_0_6_res_" + num2str(max_res) +"_dx_10_dz_10_resus_rate_0_05_width_500_depth_20_CONVSIM"");"); % non stratified, fixed resus 0.05 width 500m
end


p = temp.p;

index = 1;
for j = 1:max_res-1
    for i = 1:max_res-j
        x_vol_vec_max(index) = p.X_vol(i,j);
        z_vol_vec_max(index) = p.Z_vol(i,j);
        index = index +1;
    end
end






for i = 1:num_res % interating over the resolutions
    
    % declaring struct
    eval("p_" + num2str(res_vec(i)) + "= struct;");
    
    
    
    % loading data into struct
    if(strafified_fixed_resus)
        eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(res_vec(i)) +"_Stratified_therm_depth_2_thickness_3_profile_100_1_10_CONV_SIM"");"); % stratified simulations
    end
    
    if(strafified_fixed_resus_500m_radius)
        eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V6_alpha_1_kbg_0_6_res_" + num2str(res_vec(i)) +"_dx_10_dz_10_resus_rate_0_05_width_500_depth_20_CONVSIM"");"); % non stratified, fixed resus 0.05 width 500m
    end
    
    if(stratified_variable_resus_type_1)
        eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(res_vec(i)) +"_Stratified_therm_depth_5_thickness_3_profile_100_1_10_variable_resus_resp_fn_1_CONVSIM"");"); % stratified simulations, variable resus type 1
    end
    
%     
%     if(stratified_variable_resus_type_2)
%         eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V6_alpha_1_kbg_0_2_res_" + num2str(res_vec(i)) +"_Stratified_therm_depth_5_thickness_3_profile_100_1_10_variable_resus_resp_fn_2_CONVSIM"");"); % stratified simulations, variable resus type 2
%     end
    
    if(stratified_variable_resus_type_3_width_500)
        eval("p_" + num2str(res_vec(i)) + "= load(""2D_benthic_results_V6_alpha_1_kbg_0_6_res_" + num2str(res_vec(i)) +"_Stratified_therm_depth_5_thickness_3_profile_100_1_10_variable_resus_resp_fn_2_CONVSIM"");"); % stratified simulations, variable resus type 2
    end
    

    
    
    
    eval("  Y_t = p_" + num2str(res_vec(i)) + ".Y_t;");
    eval("    p = p_" + num2str(res_vec(i)) + ".p;");
    
    
    
    A = Y_t(end,1:sum(1:p.Xn-1)); % phytoplankton
    Rd = Y_t(end,sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
    D = Y_t(end,2*sum(1:p.Xn-1) +1 : 3*sum(1:p.Xn-1)); % ditritus
    Rs = Y_t(end,3*sum(1:p.Xn-1)+1 : 3*sum(1:p.Xn-1)+ p.Xn-1); % sedimented nutrients
    B = Y_t(end, 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae
    
    %Calculating nutrient content at t = Tend
    n_algae_end     = p.volumes_cyl.*A'*p.q;
    n_dissolved_end = p.volumes_cyl.*Rd';
    n_detritus_end  = p.volumes_cyl.*D';
    n_sediment_end  = Rs.*p.Area_bottom_cyl;
    n_benthic_end   = B.*p.q.*p.Area_bottom_cyl;
    ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);
    
    ntot_vec(i) = ntot_end;
    
    % storing total nutrient concentration for each state variable
    detritus_tot(i) = sum(sum(n_detritus_end));
    dissolved_tot(i) =  sum(sum(n_dissolved_end));
    algae_tot(i) = sum(sum(n_algae_end));
    sed_tot(i) = sum(n_sediment_end);
    benth_tot(i) = sum(n_benthic_end);
    
    %  A_interp_vec = griddata(x_vol_vec,z_vol_vec, A,x_vec,z_vec, 'nearest');
    
    
    % flattened coordinates of grid element centers
    x_vol_vec = zeros(1,length(A));
    z_vol_vec = zeros(1,length(A));
    
    index = 1;
    for j = 1:p.Xn-1
        for ii = 1:p.Zn-j
            x_vol_vec(index) = p.X_vol(ii,j);
            z_vol_vec(index) = p.Z_vol(ii,j);
            index = index +1;
        end
    end
    
    
    % Interpolating to a refined grid (same grid as the highest resolution)
    query_points = eval("linspace(0,p_" + num2str(res_vec(i)) + ".p.W,max_res-1);");
    sed_points = eval("zeros(1, p_"+num2str(res_vec(i))+".p.Xn-1);");
    sed_points(:) = eval("( p_"+num2str(res_vec(i))+".p.X(1,2:end) +  p_"+num2str(res_vec(i))+".p.X(1,1:end-1))/2;");
    
    eval("algae_matrix("+num2str(i)+",:) = griddata(x_vol_vec, z_vol_vec, A.*p_"+ num2str(res_vec(i))+".p.q,x_vol_vec_max,z_vol_vec_max, 'v4');");
    eval("detr_matrix("+num2str(i)+",:) = griddata(x_vol_vec, z_vol_vec, D ,x_vol_vec_max,z_vol_vec_max, 'v4');");
    eval("nutr_matrix("+num2str(i)+",:) = griddata(x_vol_vec, z_vol_vec, Rd,x_vol_vec_max,z_vol_vec_max, 'v4');");
    eval("sed_vec(i,:) = interp1(sed_points,Rs, query_points,'pchip');");
    eval("benth_vec(i,:) = interp1(sed_points,B, query_points,'pchip');");
end

%% %%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%
close all

% sediment nutrients
figure(1)
for i =5:length(res_vec) % i loops over each point in the sediment layer, i=1 corresponds to the leftmost point at the center of the lake.
    %plot(res_vec, abs((sed_vec(:,i)-sed_vec(end,i))./sed_vec(end,i)));
    plot(res_vec',sed_vec(:,i)/(sed_vec(end,i)));
    hold on
end
title("Local comparison of the sediment nutrient content (normalized).");
hold off

% benthic algae
figure(2)
for i =5:length(res_vec) % i loops over each point in the benthic layer,  i=1 corresponds to the leftmost point at the center of the lake.
    plot(res_vec', benth_vec(:,i)./benth_vec(end,i));
    hold on
end
title("Local comparison of the benthic nutrient content (normalized).");
hold off


% plots of local convergence detritus
figure(3)
for i =1:sum(1:max_res-1)
    plot(res_vec', detr_matrix(:,i)./detr_matrix(end,i));
    hold on
end
title("Local comparison of the detritus nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of local convergence dissolved nutrients
figure(4)
for i =1:sum(1:max_res-1)
    plot(res_vec', nutr_matrix(:,i)./nutr_matrix(end,i));
    hold on
end
title("Local comparison of the dissolved nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of local convergence dissolved nutrients
figure(5)
for i =1:sum(1:max_res-1)
    plot(res_vec, algae_matrix(:,i)./algae_matrix(end,i));
    hold on
end
title("Local comparison of the plankton nutrient content (normalized).");
xlabel("resolution");
ylabel("Concentration (normalized to highest res)");
hold off

% plots of the total nutrient content percentage of each state variable
figure(6)
plot(res_vec, detritus_tot./detritus_tot(end), res_vec, dissolved_tot./dissolved_tot(end), res_vec, algae_tot./algae_tot(end));
%plot(res_vec, detritus_tot./detritus_tot(end), res_vec, dissolved_tot./dissolved_tot(end));
hold on
plot(res_vec, benth_tot./benth_tot(end), res_vec, sed_tot./sed_tot(end));
%legend("detritus","dissolved nutrients","algae","benthic algae","sediment");
title("Normalized total nutrient of each state variable");
legend("Detritus", "dissolved nutr.","pel. algae","benth. algae","sediment");
hold off









