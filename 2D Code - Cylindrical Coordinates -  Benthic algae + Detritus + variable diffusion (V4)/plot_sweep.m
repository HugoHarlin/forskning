%% SCRIPT INFO
% Plots concentrations at steady state, either in indivudual plots or
% in one compound figure using subplot.


for kbg_index = 2
    % for resusp_index = 1:2
    for remin_index = 3
        for dx = [1, 10, 100] %[1,10,100]
            for dz = [1, 10, 100] %, [1,10,100]
                
                Xn = "16"; % horizontal resolution -1;
                Zn = "16"; % vertical resolution -1;
                alpha_str = "1_5";
                kbg_str = ["0","0_2","0_4","0_8","2_0"];
                resusp_str = ["0_0", "0_1"];
                remin_str = ["0_01","0_05","0_1"];
                % filename = "2D_benthic_results_V4_dx_" + num2str(dx) + "_dz_" + num2str(dz) + "_Xn_" + Xn + "_Zn_" + Zn + "_alpha_" + alpha_str + "_kgb_"+ kbg_str(kbg_index);
                % filename = "2D_benthic_results_V4_dx_" + num2str(dx) + "_dz_" + num2str(dz) + "_Xn_" + Xn + "_Zn_" + Zn + "_alpha_" + alpha_str + "_kgb_"+ kbg_str(kbg_index) + "_resuspRate_" +resusp_str(resusp_index);
                filename = "2D_benthic_results_V4_dx_" + num2str(dx) + "_dz_" + num2str(dz) + "_Xn_" + Xn + "_Zn_" + Zn + "_alpha_" + alpha_str + "_kbg_"+ kbg_str(kbg_index) +  "_reminrate_"+ remin_str(remin_index) ;
                
                load(filename);
                plot_results_fn(t,Y_t,p);
            end
        end
    end
end

