% Automatic plotting script.

% This script loads the data files (.mat) in the current folder and runs the plot script.
clc
clear 
close all

% we add the path of the autoplot script and the plot_results_v6 script
addpath('C:\Users\huha0015\Documents\PhD\forskning\2D Code - Cylindrical Coordinates -  Benthic algae + Detritus + variable diffusion (V6) - square grid - phytoplankton = benthic algae\autoplot_script.m');
addpath("C:\Users\huha0015\Documents\PhD\forskning\2D Code - Cylindrical Coordinates -  Benthic algae + Detritus + variable diffusion (V6) - square grid - phytoplankton = benthic algae\plot_results_fn_v6.m")


% we store all the files in the current directory.
all_files = dir();


%% running standard plot script
if(false)
    % looping over each file in the current directory
    for i=1:length(all_files)
        
        name_temp = all_files(i).name; % name of file i
        if(length(name_temp) >4) % the file has to be at least 5 characters long to be a .mat file
            
            if(strcmp(name_temp(end-3:end),'.mat')) % the file is a .mat file, and we want to load it into the workspace.
                
                load(name_temp); % we load the .mat file
                
                if(isfield(p,'model_version') == 0)
                    p.model_version = 6; % if the model version parameter doens't exist (as in older sims) we hardcode it. Make sure it is correct!
                end
                
                plot_results_fn_v6(t,Y_t,p); % we plot the results, which are saved in the working directory.
                clear t Y_t p % we remove the current data from the workspace before loading the next.
                
            end
            
            
        end
    end
end

%% saving state variables 
    % looping over each file in the current directory
    res_tot = struct([]); % array of structs, where each struct containts the extracted results from each index
    index = 1; % index for keeping track if the number of simulation files loaded and stored
    for i=1:length(all_files)
        
        name_temp = all_files(i).name; % name of file i
        if(length(name_temp) >4) % the file has to be at least 5 characters long to be a .mat file
            
            if(strcmp(name_temp(end-3:end),'.mat')) % the file is a .mat file, and we want to load it into the workspace.
                
                load(name_temp); % we load the .mat file
                
                if(isfield(p,'model_version') == 0)
                    p.model_version = 6; % if the model version parameter doens't exist (as in older sims) we hardcode it. Make sure it is correct!
                end
                
                
                % we save the results from each run in the vector
                [r, A,Rd,D,Rs,B,ntot_0,n_algae_end,n_dissolved_end,n_detritus_end,n_sediment_end,n_benthic_end,ntot_end] = extract_results_v6(t,Y_t,p);
                
                
                % saving results in struct vecotor res_tot.
                
                % state variables at steady state
                res_tot(index).A  = A;
                res_tot(index).Rd = Rd;
                res_tot(index).D  = D;
                res_tot(index).Rs = Rs;
                res_tot(index).B  = B;
                
                % nutrients at steady state
                res_tot(index).ntot_0   = ntot_0;
                res_tot(index).ntot_end = ntot_end;
                
                res_tot(index).n_algae_end     = n_algae_end;
                res_tot(index).n_dissolved_end = n_dissolved_end;
                res_tot(index).n_detritus_end  = n_detritus_end;
                res_tot(index).n_sediment_end  = n_sediment_end;
                res_tot(index).n_benthic_end   = n_benthic_end;
                
                % resuspension coeff
                 res_tot(index).r = r;
                
                
                index = index +1;
               % plot_results_fn_v6(t,Y_t,p); % we plot the results, which are saved in the working directory.
                clear r
                clear t Y_t p % we remove the current data from the workspace before loading the next.
                clear A Rd D Rs B 
                clear ntot_0 ntot_end
                clear n_algae_end n_dissolved_end n_detritus_end n_sediment_end n_benthic_end 
            end
            
            
        end
    end

%% plotting lakewide nutrient portions as a function of resupension coefficient
num_sims = length(res_tot(:));
resus_vec = zeros(1,num_sims);
tot_alg = zeros(1,num_sims); % percentage of nutrients in algae
tot_benth = zeros(1,num_sims); % percentage of nutrients in benthic algae

for i=1:num_sims
    resus_vec(i) = res_tot(i).r;
    tot_alg(i) = sum(res_tot(i).n_algae_end)/res_tot(i).ntot_end;
    tot_benth(i) = sum(res_tot(i).n_benthic_end)/res_tot(i).ntot_end;
end

%% plots of nutrient portions as a function of resuspension

close all
plot(resus_vec, tot_benth+tot_alg,'lineWidth',2, 'Color', 'black')
hold on
%resus_vec = sort(resus_vec);
%tot_benth = sort(tot_benth);
plot(resus_vec, tot_benth,'lineWidth',2)
hold on
%tot_alg = sort(tot_alg);
plot(resus_vec, tot_alg,'lineWidth',2)

xlabel('Resuspension coefficient')
ylabel('Portion of total nutrients')
legend('Total biomass','Pelagic algae','Benthic algae', 'Box','off');
%ylim([0 0.6]);
set(gca,'FontSize',22)
set(gca,'Box','off');
