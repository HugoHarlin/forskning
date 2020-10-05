
function [] = plot_results_fn_v6(t,Y_t,p)
% plots simulation results in one figure and saves them as a jpg.

%% Extracting results
ntot_0 = p.ntot_0;


% Time of plotted results
%time = length(t); % winter
%time = length(t)-183; % summer
%time = length(t)-274; % spring
%time = length(t)-91; % fall
time = length(t);

A = Y_t(time,1:sum(1:p.Xn-1)); % phytoplankton
Rd = Y_t(time,sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
D = Y_t(time,2*sum(1:p.Xn-1) +1 : 3*sum(1:p.Xn-1)); % ditritus
Rs = Y_t(time,3*sum(1:p.Xn-1)+1 : 3*sum(1:p.Xn-1)+ p.Xn-1); % sedimented nutrients
B = Y_t(time, 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae


%Calculating nutrient content at t = Tend
n_algae_end     = p.volumes_cyl.*A'*p.q;
n_dissolved_end = p.volumes_cyl.*Rd';
n_detritus_end  = p.volumes_cyl.*D';
n_sediment_end  = Rs.*p.Area_bottom_cyl;
n_benthic_end   = B.*p.q.*p.Area_bottom_cyl;
ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);

% creation of refined grid for interpolation.
%[X_vol_new,Y_vol_new] = grid_interpolation_fn(p,100,100);


%% Calculation of the light intensity
I = zeros(sum(1:p.Xn-1),1);

index = 1;
temp = 1;
if(isfield(p,'seasonality')) % checks if the 'seasonality' variable exists
    if(p.seasonality)
        if(p.seasonality_light) % if seasonal light is active we need to change p.I0 to the correct value
            % seasonal light is modelled with a sin-function that varies from
            % p.minLight to p.maxlight
            p.I0 = (p.maxLight-p.minLight)/2*sin(2*pi*(90+mod(t(time),365))/365) + (p.minLight + p.maxLight)/2;
        end
        
        if(isfield(p,'seasonality_thermoC'))
            if(p.seasonality_thermoC) % varying depth of the thermocline
                p.thermocline_depth = (p.maxTherm -p.minTherm)/2*sin(2*pi*(90+mod(t(time),365))/365) + (p.minTherm + p.maxTherm )/2;
            end
        end
        
    end
end

for j = 1:p.Xn-1 % pooling over each column (j = column nr)
    int = 0;
    for i = 1:p.Zn - temp % traversing down one column (i = row nr)
        
        dz_middle = p.Z_vol(i,j) - p.Z(i,j); % distance from bottom of grid element (i-1,j) to the center of grid element (i,j)
        
        if(i==1) % width of grid element (i,j) in the z-direction
            dz = p.Z(i+1,j);
        else
            dz = p.Z(i+1,j) - p.Z(i,j);
        end
        
        int_middle = int + (p.kD.*D(index) +p.kA*A(index)+ p.kbg)*dz_middle; % this is the integral down to the center of grid element (i,j).
        int = int + (p.kD.*D(index) +p.kA*A(index)+ p.kbg)*dz; % here we integrate to the bottom of grid element (i,j)
        I(index) =  p.I0 .* exp( - int_middle);
        index = index +1;
    end
    temp = temp +1;
end
p.I = I;

%% Reformatting state variables and light into matrixes
I_matrix  = NaN.*ones(p.Zn-1, p.Xn-1);
A_matrix  = NaN.*ones(p.Zn-1, p.Xn-1);
Rd_matrix = NaN.*ones(p.Zn-1, p.Xn-1);
D_matrix  = NaN.*ones(p.Zn-1, p.Xn-1);

temp = 1;
index = 1;
for j = 1:p.Xn-1 % looping over each column
    for i = 1:p.Zn - temp % traversing down one column (i = row nr)
        I_matrix(i,j) = I(index);
        A_matrix(i,j)  = A(index);
        Rd_matrix(i,j) = Rd(index);
        D_matrix(i,j)  = D(index);
        index = index +1;
    end
    temp = temp +1;
end

% flattened coordinates of grid element centers
x_vol_vec = zeros(1,length(A));
z_vol_vec = zeros(1,length(A));

index = 1;
for j = 1:p.Xn-1
    for i = 1:p.Zn-j
        x_vol_vec(index) = p.X_vol(i,j);
        z_vol_vec(index) = p.Z_vol(i,j);
        index = index +1;
    end
end


% flattened coordinates of grid element centers
%x_vec = zeros(1,sum(2:p.Xn)+p.Xn);
%z_vec = zeros(1,sum(2:p.Xn)+p.Xn);


% X_vec and Z_vec are the coordinates of the grid corneres (nodes),
% flattened in a matrix column by column, in decending order starting from
% the surface.

x_vec(1:p.Xn) =  p.X(:,1);
z_vec(1:p.Xn) =  p.Z(:,1);
index = p.Xn+1;
for j = 2:p.Xn
    for i = 1:p.Zn-j+2
        x_vec(index) = p.X(i,j);
        z_vec(index) = p.Z(i,j);
        index = index +1;
    end
end

%% interpolating bulk state variables for plotting

A_interp_vec = griddata(x_vol_vec,z_vol_vec, A,x_vec,z_vec, 'v4'); % nearest
Rd_interp_vec = griddata(x_vol_vec,z_vol_vec, Rd,x_vec,z_vec, 'v4'); % nearest
D_interp_vec = griddata(x_vol_vec,z_vol_vec, D,x_vec,z_vec, 'v4'); % nearest
A_interp = NaN*ones(p.Xn);
Rd_interp = NaN*ones(p.Xn);
D_interp = NaN*ones(p.Xn);



%  A_interp(:,1) = A_interp_vec(1:p.Xn);
% Rd_interp(:,1) = Rd_interp_vec(1:p.Xn);
% D_interp(:,1) = D_interp_vec(1:p.Xn);
% index = p.Xn+1;

j = 1:p.Xn;
i = 1;
index = 1:p.Xn;

A_interp(j,i) = A_interp_vec(index);
Rd_interp(j,i) = Rd_interp_vec(index);
D_interp(j,i) = D_interp_vec(index);


index =p.Xn+1; % p.Xn+1;
for i =2:p.Xn
    for j=1:p.Zn-i+2
        A_interp(j,i) = A_interp_vec(index);
        Rd_interp(j,i) = Rd_interp_vec(index);
        D_interp(j,i) = D_interp_vec(index);
        index = index +1;
    end
end
test = 1;

%% Figure index
figureIndex = 1;
close all

%% Plot of simulation results

if(true)
    font_size = 18;
    
    %% plots of all concentrations in one figure using griddata
    if(true)
        figure(figureIndex)
        figureIndex = figureIndex +1;
        %  set(gcf, 'Position',  [1700, -660,  1730, 1300]) % for widescreen
        set(gcf, 'Position',  [1700, 50,  1730, 1300]) % for widescreen
        % set(gcf, 'Position',  [0, 0, 1200, 1000]) % if on laptop
        tile_fig =  tiledlayout(3,3,'TileSpacing','Compact');
        tile_fig.Padding = 'compact';
        title_str = "";
        %title_str = "dx: " + num2str(p.dx(1,1)) + "  dz: " + num2str(p.dz(1,1)) + "  Xn: " + num2str(p.Xn) + "  Zn: " + num2str(p.Zn) + "  alpha: " + num2str(p.alpha) + "  kbg: " + num2str(p.kbg) + "";
        
        if(p.stratified)
            title_str = "Stratified. depth: " + num2str(p.thermocline_depth) + ", th: " + num2str(p.thermocline_thickness) + ", " ;
        end
        
        title_str = title_str + "resolution: " + num2str(p.Xn-1) + " kbg: " + num2str(p.kbg) + ""; %  + " alpha: " + num2str(p.alpha)
        
        if(p.constant_resuspension == 1)
            title_str = title_str + ", resus rate: " + num2str(p.resus(1)) ;
        else
            if(p.response_type == 1 ||  p.response_type == 2 ||  p.response_type == 3)
                title_str = title_str +", var. resus type " + num2str(p.response_type);
            end
            
            if(p.stepFun_resus)
                title_str = title_str +", stepfun resus, " + num2str(p.upper_thresh) + ", " + num2str(p.lower_thresh) + ", " + num2str(p.upper_resusp) + ", " + num2str(p.  lower_resusp);
            end
            
            
            if(p.manual_in_therm_resus)
                title_str = title_str +  " man. resus " + num2str(p.manual_therm_resus_val);
            end
        end
        
        if(~p.stratified) % not stratified
            title_str = title_str + " dx: " + num2str(p.dx(1,1));
            title_str = title_str + " dz: " + num2str(p.dz(1,1));
        end
        
        title(tile_fig,title_str,'FontSize', font_size + 12); % shared title for all plots
        x_aspect = 1.2;
        
        
        % extracting diagonal elements for triangle elements
        
        A_diag = zeros(3*(p.Xn-1),1);
        Rd_diag = zeros(3*(p.Xn-1),1);
        D_diag = zeros(3*(p.Xn-1),1);
        I_diag = zeros(3*(p.Xn-1),1);
        
        i = p.Zn-1;
        index = 1;
        
        for j = 1: p.Xn-1
            A_diag(index:index+2)  = A_matrix(i,j);
            Rd_diag(index:index+2) = Rd_matrix(i,j);
            D_diag(index:index+2)  = D_matrix(i,j);
            I_diag(index:index+2)  =  I_matrix(i,j);
            
            i = i -1;
            index = index +3;
        end
        
        
        % connectivity matrix for the triangle elements of the diagonal.
        T = zeros(p.Xn-1,3);
        counter = 1;
        for i = 1:p.Xn-1
            T(i,:) = 1;
            T(i,1) = counter;
            T(i,2) = counter+1;
            T(i,3) = counter+2;
            counter = counter +3;
        end
        
        % x and z coordinates of the points of the triangle elements of the
        % diagonal
        x_triag = zeros(3*(p.Xn-1),1); % for trisurf to work the vector containing the x-coordinates must be of sufficient length to store unique valeus for all nodes, regardless of values are shared or not.
        z_triag = zeros(3*(p.Xn-1),1);
        
        i = p.Zn;
        idx = 1;
        for j = 1: p.Xn-1
            x_triag(idx) = p.X(i,j);
            x_triag(idx+1) = p.X(i,j);
            x_triag(idx+2) = p.X(i,j+1);
            
            z_triag(idx) = p.Z(i,j);
            z_triag(idx+1) = p.Z(i-1,j);
            z_triag(idx+2) = p.Z(i-1,j);
            
            idx = idx +3;
            i = i -1;
        end
        
        
        test = 1;
        
        %% phytoplankton
        
        % Plots algal density where the values have been interpolated to
        % produce a smooth figure.
        nexttile
        % plotting triangular elements
        surf(p.X,p.Z,p.q.*A_interp, 'edgecolor','none') % plotting interpolated values of A
        hold on
        
        set(gca, 'Ydir', 'reverse')
        caxis([0 inf]);
        shading(gca,'interp')
        grid off
        az =0;
        el = 90;
        view(az,el);
        hold on
        
        
        % plotting white triangle to cover diagonal sawtooth
        trisurf([1 2 3], [0 p.W p.W], [p.Lmax p.Lmax 0],(1.01*max(max(p.q.*A_interp))).*[1 1 1], 'edgecolor','none', 'FaceColor','white')
        hold on
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        if(p.stratified)
            value = 1.01*max(max(p.q.*A_interp));
            therm =  surf([0 0; p.W p.W],[p.thermocline_depth p.thermocline_depth + p.thermocline_thickness; p.thermocline_depth p.thermocline_depth + p.thermocline_thickness],[value value; value value]);
            therm.EdgeColor = 'none';
            therm.FaceAlpha  = 0.1;
            therm.FaceColor = 'black';
        end
        
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
        if(z_theory_algae> 0 && z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[max(max(p.q.*A))+1, max(max(p.q.*A))+1], '--','color','black','LineWidth',2); % plotting the max depth line
            mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            h = annotation('textbox','String',mystr,'fontsize',font_size -2 ,'EdgeColor','none');
            h.Position=[.115 .635 .3581 .0879];
        end
        
        hold off
        
        %% Dissolved nutrients
        
        nexttile
        
        %surf(p.X_vol, p.Z_vol, Rd_matrix,  'edgecolor','none');
        
        surf(p.X,p.Z,Rd_interp, 'edgecolor','none') % plotting interpolated values of Rd
        hold on
        
        
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
        hold on
        
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[1.01*max(max(A)), 1.01*max(max(A))], '--','color','black','LineWidth',2); % plotting the max depth line
        end
        hold on
        % plotting a white triangle covering the portruding sawtooth along
        % the diagonal.
        trisurf([1 2 3], [0 p.W p.W], [p.Lmax p.Lmax 0],(1.01*max(max(Rd_interp))).*[1 1 1], 'edgecolor','none', 'FaceColor','white') % plotting the diagonal elements
        hold on
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        if(p.stratified)
            value = 1.01*max(max(Rd_interp));
            therm =  surf([0 0; p.W p.W],[p.thermocline_depth p.thermocline_depth + p.thermocline_thickness; p.thermocline_depth p.thermocline_depth + p.thermocline_thickness], [value value; value value]);
            therm.EdgeColor = 'none';
            therm.FaceAlpha  = 0.1;
            therm.FaceColor = 'black';
        end
        
        %ylabel("Depth [m]");
        %xlabel("distance from lake center [m]");
        hold off
        
        %% Detritus
        nexttile
        %surf(p.X_vol, p.Z_vol, D_matrix,  'edgecolor','none');
        
        surf(p.X,p.Z,D_interp, 'edgecolor','none') % plotting interpolated values of D
        hold on
        
        
        % plots the max survival depth of pelagic algae as a dotted line
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[max(max(A))+1, max(max(A))+1], '--','color','black','LineWidth',2); % plotting the max depth line
        end
        
        shading interp
        
        % plots a clear triangular patch covering the sawtooth diagonal
        trisurf([1 2 3], [0 p.W p.W], [p.Lmax p.Lmax 0],(1.01*max(max(D_interp))).*[1 1 1], 'edgecolor','none', 'FaceColor','white') % plotting the diagonal elements
        hold on
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        if(p.stratified)
            value = 1.01*max(max(D_interp));
            therm =  surf([0 0; p.W p.W],[p.thermocline_depth p.thermocline_depth + p.thermocline_thickness; p.thermocline_depth p.thermocline_depth + p.thermocline_thickness], [value value; value value]);
            therm.EdgeColor = 'none';
            therm.FaceAlpha  = 0.1;
            therm.FaceColor = 'black';
        end
        
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
        hold off
        
        %% Sediment nutrients
        sed_points = zeros(1,p.Xn-1); % creating points for sediment plot
        sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        query_points = linspace(0,p.W,200);
        %sed_interp = interp1(sed_points,Rs, query_points,'makima');
        sed_interp = interp1(sed_points,Rs, query_points,'pchip');
        nexttile
        plot(query_points,sed_interp,'color','black');
        
        
        max_survival_depth = -1./p.kbg .* log(p.H_benth./( p.I0.*(p.Gmax_benth./p.lbg_benth -1) )); % maximum survival depth of benthic algae assuming no phytoplankton
        max_survival_depth_x = p.W*((max_survival_depth./p.Lmax -1)*(p.Lmax/(p.Lmin-p.Lmax)))^(1/p.alpha); % corressponding distance from the lake center
        if(max_survival_depth_x > 0)
            xline(max_survival_depth_x,'LineStyle','--','LineWidth',2) % maximum survival depth
        end
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        if(p.stratified)
            x1 = (p.W/p.Lmax)*(p.Lmax - p.thermocline_depth - p.thermocline_thickness);
            x2 = (p.W/p.Lmax)*(p.Lmax - p.thermocline_depth);
            axes = gca;
            y =max(axes.YLim);
            hold on
            temp = area([x1, x2],[y y]);
            
            % the max y values is increased when the bar representing the
            % thermocline is added (this is done by default so the line isn't
            % right at the top edge of the plot. But, ths is exactly what we want,
            % so we change the value back to before we added the shaded area for the thermocline.
            axes.YLim(end) = y;
            temp.EdgeColor = 'none';
            temp.FaceAlpha  = 0.1;
            temp.FaceColor = 'black';
        end
        
        
        
        pbaspect([x_aspect 1 1])
        title("Sediment Nutrients [mgP/m^2]");
        ylabel("area density [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        set(gca,'fontSize',font_size);
        
        %% Benthic algae
        % creating points for benthic algae plot
        benth_points = zeros(1,p.Xn-1);
        benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
        new_res = 300; % resolution of interpolated plot
        query_points = linspace(0,p.W,new_res);
        %benth_interp = interp1(sed_points,B, query_points,'makima');
        benth_interp = interp1(sed_points,B.*p.q, query_points,'pchip','extrap'); % interpolated net growtht values, unit [mgP/(m^2 day)]
        fig = nexttile;
        
        colororder({'b','black'});
        
        % calculating the light/nutrient limitation for the benthic algae
        benthic_limitation = zeros(1,new_res);
        
        
        % extracting dissolved nutrient concentrations and light intensities along the bottom
        Rd_bottom_temp = zeros(1,p.Xn-1);
        I_bottom_temp = zeros(1,p.Xn-1);
        benth_prod_interp = zeros(1,new_res);
        
        
        b_idx = p.Zn-1;
        for i=1:p.Xn-1
            Rd_bottom_temp(i) = Rd(b_idx);
            I_bottom_temp(i) = I(b_idx);
            b_idx = b_idx + p.Zn-1 -i;
        end
        
        % interpolating to the refined grid used in the plot
        % benthic_limitation_interp =   interp1(sed_points, benthic_limitation, query_points,'pchip','extrap');
        % benth_prod_interp = interp1(sed_points, benth_prod, query_points,'pchip','extrap');
        Benth_interp = interp1(sed_points, B, query_points,'pchip','extrap'); % spline
        %Benth_interp = interp1(sed_points, B, query_points,'spline','extrap'); % spline
        Rd_bottom_interp = interp1(sed_points, Rd_bottom_temp, query_points,'pchip','extrap');
        I_bottom_interp = interp1(sed_points,  I_bottom_temp, query_points,'pchip','extrap');
        
        nutrient_limited_growth = zeros(1,new_res);
        light_limited_growth = zeros(1,new_res);
        
        for i=1:new_res
            nutrient_limited_growth(i) =  p.Gmax_benth.*Benth_interp(i).*(Rd_bottom_interp(i)./(Rd_bottom_interp(i)+p.M_benth));
            light_limited_growth(i) = p.Gmax_benth./p.kB.*log((p.H_benth + I_bottom_interp(i))./(p.H_benth + I_bottom_interp(i).*exp(-p.kB.*Benth_interp(i))));
            
            % benth_prod(i) = min(nutrient_limited_growth,light_limited_growth) -(p.lbg_benth + p.resus_benth)*B(i);
            
            % Due to interpolation the growths can be negative in extreme
            % cases, of so we set the growth to zero in order not to mess
            % upp the plots.
            if((nutrient_limited_growth(i) < 0)) nutrient_limited_growth(i) = 0; end
            if((light_limited_growth(i) < 0)) light_limited_growth(i) = 0; end
            
            
            
            %   benth_prod_interp(i) = min(nutrient_limited_growth(i),light_limited_growth(i))/Benth_interp(i);
            
            
            
            if(Benth_interp(i) > 1e-8) % threshold value, below which the benthic algae is considered extinct and we don't compute the growth limitation (which requires algae present. if not done the specific growth plot looks very erratic where there is miniscule/no ammounts of benthic algae.
                benth_prod_interp(i) = min(nutrient_limited_growth(i),light_limited_growth(i))/Benth_interp(i);
            else
                benth_prod_interp(i) = 0; % the calculation of the benthic growth limitation requires there to be benthic algae present since the benthic algae shadows itself and this is taken into account.
                % therefore the benthic production is set to zero when the
                % benthic algae is locally extinct so that the observer
                % isn't confused by the erratic specific production where
                % there is practically zero benthic algae.
            end
            
            
            if(nutrient_limited_growth(i) > light_limited_growth(i))
                
                benthic_limitation(i) = 1; % benthic algae is light limited
            else
                benthic_limitation(i) = -1; % benthic algae is nutrient limited
            end
        end
        benthic_limitation_interp =  benthic_limitation;
        
        
        % old code, where the nutrient/light limitation switch was
        % calculated before interpolating.
        %         for i=1:p.Xn-1
        %             nutrient_limited_growth =  p.Gmax_benth.*B(i)'.*(Rd(b_idx)./(Rd(b_idx)+p.M_benth));
        %             light_limited_growth = p.Gmax_benth./p.kB.*log((p.H_benth + I(b_idx))./(p.H_benth + I(b_idx).*exp(-p.kB.*B(i))));
        %
        %             % benth_prod(i) = min(nutrient_limited_growth,light_limited_growth) -(p.lbg_benth + p.resus_benth)*B(i);
        %             benth_prod(i) = min(nutrient_limited_growth,light_limited_growth)/B(i);
        %
        %             if(nutrient_limited_growth > light_limited_growth)
        %
        %                 benthic_limitation(i) = 1; % benthic algae is light limited
        %             else
        %                 benthic_limitation(i) = -1; % benthic algae is nutrient limited
        %             end
        %
        %             b_idx = b_idx + p.Zn-1 -i;
        %         end
        
        
        
        lim_switch_index = 0;
        index = 1;
        % not the prettiest but get the job done...
        for i=1:length(benthic_limitation_interp)-1
            if(  sign(benthic_limitation_interp(i)) ~= sign(benthic_limitation_interp(i+1)) )
                lim_switch_index(index) = i;
                index = index +1;
            end
        end
        
        
        
        
        if(lim_switch_index ~= 0)
            color = 'bl';
            start_index = 1;
            for i=1:length(lim_switch_index)
                if(sign(benthic_limitation_interp(lim_switch_index(i)-1)) < 0)
                    
                    % plotting the benthic algae that is light limited [mgP/m^3]
                    color = 'red';
                else
                    
                    % plotting the benthic algae that is nutrient limited
                    color = 'blue';
                end
                hold on
                yyaxis left
                plot(query_points(start_index:lim_switch_index(i)), benth_interp(start_index:lim_switch_index(i)),color, 'LineWidth', 2,'LineStyle','-','Marker', 'none');
                start_index = lim_switch_index(i);
            end
            hold on
            if(sign(benthic_limitation_interp(lim_switch_index(i)+1)) < 0)
                
                % plotting the benthic algae that is light limited [mgP/m^3]
                color = 'red';
            else
                % plotting the benthic algae that is nutrient limited
                color = 'blue';
            end
            yyaxis left
            plot(query_points(lim_switch_index(end):end), benth_interp(lim_switch_index(end):end),'-', 'color', color, 'LineWidth', 2);
        else
            if(benthic_limitation_interp(end) == -1)
                color = 'red';
            else
                color = 'blue';
            end
            yyaxis left
            plot(query_points, benth_interp,color,'LineWidth', 2,'LineStyle','-');
        end
        
        % plotting the maximum survival depth
        if(max_survival_depth_x > 0)
            xline(max_survival_depth_x,'LineStyle','--','LineWidth',2)
        end
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        if(p.stratified)
            x1 = (p.W/p.Lmax)*(p.Lmax - p.thermocline_depth - p.thermocline_thickness);
            x2 = (p.W/p.Lmax)*(p.Lmax - p.thermocline_depth);
            axes = gca;
            y =max(axes.YLim);
            hold on
            temp = area([x1, x2],[y y]);
            
            % the max y values is increased when the bar representing the
            % thermocline is added (this is done by default so the line isn't
            % right at the top edge of the plot. But, ths is exactly what we want,
            % so we change the value back to before we added the shaded area for the thermocline.
            axes.YLim(end) = y;
            temp.EdgeColor = 'none';
            temp.FaceAlpha  = 0.1;
            temp.FaceColor = 'black';
        end
        
        
        % pbaspect([x_aspect 1 1])
        %title("Benthic Algae [mgP/m^2]");
        ylabel("area density [mgP/m^2]");
        xlabel("distance from lake center [m]");
        ylim([0 inf]);
        
        
        % we add a second plot of the specific production on top of the
        % benthic nutrient content, adding a second y-axis to the right
        hold on
        
        % right axis
        yyaxis right
        plot(query_points, benth_prod_interp,'color','black');
        ylabel("benthic specific production");
         ylabel("specific production");
        
        set(gca,'fontSize',font_size);
        pbaspect([x_aspect 1 1])
        
        %% light/nutrient limitation
        % calculation  of light/nutrient limitation at each point in lake.
        
        new_res = 200;
        X_new = zeros(new_res+1); % Mesh-spacing in x-dimension
        Z_new = zeros(new_res+1); % Mesh-spacing in y-dimension
        
        for ii=1:new_res+1
            X_new(ii,:) = (0:p.W/(new_res):p.W)'; % even spacing horizontally
        end
        
        for ii=1:new_res+1
            Z_new(:,ii) = (0:p.Lmax/(new_res):p.Lmax)';  % even spacing vertically
        end
        
        
        % coordinates of the center of each grid element
        X_vol_new = zeros(new_res);
        Z_vol_new = zeros(new_res);
        
        for jjj=1:new_res
            for ii =1:new_res
                x_idx = jjj+1; % Required (static): When using the nested for-loop variable for indexing a sliced array, you must use the variable in plain form, not as part of an expression.
                z_idx = ii+1;
                X_vol_new(ii,jjj) = (X_new(ii,jjj) + X_new(z_idx,jjj) + X_new(ii,x_idx) + X_new(z_idx,x_idx) ) /4;
                Z_vol_new(ii,jjj) = (Z_new(ii,jjj) + Z_new(z_idx,jjj) + Z_new(ii,x_idx) + Z_new(z_idx,x_idx) ) /4;
            end
        end
        
        
        % removing nodes outside the grid.
        jjjj = new_res;
        for ii=2:new_res
            Z_vol_new(ii,jjjj:end) = NaN;
            X_vol_new(ii,jjjj:end) = NaN;
            jjjj = jjjj -1;
        end
        
        
        % Calculating light and nutrient functional responses
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        
        l_lim_interp = griddata(p.X_vol(~isnan(p.X_vol)),p.Z_vol(~isnan(p.Z_vol)), l_lim , X_vol_new(~isnan(X_vol_new)),Z_vol_new(~isnan(Z_vol_new)), 'v4'); % interpolates to refined grid.
        n_lim_interp = griddata(p.X_vol(~isnan(p.X_vol)),p.Z_vol(~isnan(p.Z_vol)), n_lim , X_vol_new(~isnan(X_vol_new)),Z_vol_new(~isnan(Z_vol_new)), 'v4'); % interpolates to refined grid.
        
        growth_vec = p.Gmax.* min(l_lim,n_lim');  % specific growth of phytoplankton [%/day]
        %growth_vec = p.Gmax.* min(l_lim_interp,n_lim_interp');  % specific growth of phytoplankton [%/day]
        phytoplankton_primary_prod =growth_vec.*A';
        
        growth_vec = griddata(p.X_vol(~isnan(p.X_vol)),p.Z_vol(~isnan(p.Z_vol)),growth_vec , X_vol_new(~isnan(X_vol_new)),Z_vol_new(~isnan(Z_vol_new)), 'v4'); % interpolates to refined grid.
        
        limitation_vec = zeros(length(l_lim_interp),1);
        
        
        val = max(max(growth_vec)); % temp value used for the limitation boundary
        % calculating if nutrient or light limited
        for i=1:length( l_lim_interp)
            if( l_lim_interp(i) >  n_lim_interp(i))
                limitation_vec(i) = val; % nutrient limited
            else
                limitation_vec(i) = val*1.0001; % light limited
            end
        end
        
        % allocating lmitation matrix
        limitation_mtrx = NaN*zeros(new_res);
        
        %allocating specific growth matrix
        growth_mtrx = NaN*zeros(new_res);
        
        % reformating limitation data into matrix for plotting
        % first column
        limitation_mtrx(:,1) = limitation_vec(1:new_res);
        growth_mtrx(:,1) = growth_vec(1:new_res);
        
        
        index = new_res;
        for j = 2:1:new_res+1-1
            for i = 1:new_res+1 -j
                limitation_mtrx(i,j) = limitation_vec(index);
                growth_mtrx(i,j) = growth_vec(index);
                index = index +1;
            end
            
        end
        
        
        
        z_theory_algae = 1/p.kbg *log((p.I0*(p.Gmax/(p.lbg_A+p.Ad) -1))/p.H); % theoretical maximum depth for benthic growth if light limited,
        % at this depth the light limited growth is equal to respiration losses and death rate.
        z_theory_algae = round( z_theory_algae,2);
        
        nexttile
        surf(X_vol_new,Z_vol_new,growth_mtrx ,  'edgecolor','none') % plotting specific growth of pelagic algae
        shading(gca,'interp')
        hold on
        
        
        grid off
        pbaspect([x_aspect 1 1])
        az =0;
        el = 90;
        view(az,el);
        set(gca, 'Ydir', 'reverse')
        colormap(jet)
        colorbar
        title("Specific Growth of Pelagic Algae");
        xlabel("distance from lake center [m]");
        ylabel("Depth [m]");
        set(gca,'fontSize',font_size-2);
        
        %         % plotting white triangle to cover diagonal sawtooth
        %         trisurf([1 2 3], [0 p.W p.W], [p.Lmax p.Lmax 0],(1.01*max(max(growth_mtrx))).*[1 1 1], 'edgecolor','none', 'FaceColor','white')
        %         hold on
        
        % plotting the border signifying the switch from light to nutrient
        % limitation as a white line
        hold on
        contour3(X_vol_new,Z_vol_new,limitation_mtrx,'color','white','linewidth',1)
        
        % if stratified, the thermocline is plotted as a gray horizontal bar
        hold on
        if(p.stratified)
            value = 1.01*max(max(growth_mtrx));
            therm =  surf([0 0; p.W p.W],[p.thermocline_depth p.thermocline_depth + p.thermocline_thickness; p.thermocline_depth p.thermocline_depth + p.thermocline_thickness], [value value; value value]);
            therm.EdgeColor = 'none';
            therm.FaceAlpha  = 0.1;
            therm.FaceColor = 'black';
        end
        
        if(z_theory_algae < p.Lmax)
            hold on
            plot3([0; p.W],[ z_theory_algae ; z_theory_algae],[2 2], '--','color','black','LineWidth',2);
            %mystr = "Max survival depth: " + num2str(z_theory_algae) + "m";
            %h = annotation('textbox','String',mystr,'fontsize',15,'EdgeColor','none');
            %h.Position=[.54 .37 .3581 .0879];
            hold off
            set(gca,'fontSize',font_size-2);
        end
        
        %% parameters and response variables
        tt =  nexttile;
        delete(tt);
        
        perc_algae = sum(sum(n_algae_end))/ntot_end;
        perc_dissolved = sum(sum(n_dissolved_end))/ntot_end;
        perc_detritus = sum(sum(n_detritus_end))/ntot_end;
        perc_sediment = sum(n_sediment_end)/ntot_end;
        perc_benthic = sum(n_benthic_end)/ntot_end;
        [x,isterm,dir] = eventfun_V6(t(time),Y_t(time,:)',p);
        
        
        
        % Calculating light and nutrient functional responses for the
        % pelagic algae
        l_lim = I./(I+p.H);
        n_lim = Rd./(Rd+p.M);
        
        
        % this vector at index i is set to 1 if the pelagic algae in cell i is
        % light limited.
        l_lim_vec_pel = zeros(1,length(l_lim));
        
        for i=1:length(l_lim_vec_pel)
            if(l_lim(i) <  n_lim(i))
                l_lim_vec_pel(i) = 1;
            end
        end
        
        
        % calculating the light/nutrient limitation for the benthic algae
        b_idx = p.Zn-1;
        
        benth_prod = zeros(p.Xn-1,1);
        benthic_limitation = zeros(p.Xn-1,1);
        
        for i=1:p.Xn-1
            nutrient_limited_growth =  p.Gmax_benth.*B(i)'.*(Rd(b_idx)./(Rd(b_idx)+p.M_benth));
            light_limited_growth = p.Gmax_benth./p.kB.*log((p.H_benth + I(b_idx))./(p.H_benth + I(b_idx).*exp(-p.kB.*B(i))));
            
            benth_prod(i) = min(nutrient_limited_growth,light_limited_growth);
            
            if(nutrient_limited_growth > light_limited_growth)
                
                benthic_limitation(i) = 1; % benthic algae is light limited
            else
                benthic_limitation(i) = 0; % benthic algae is nutrient limited
            end
            
            b_idx = b_idx + p.Zn-1 -i;
        end
        
        
        % calculating the flux of pelagic algae into the benthos (sinking
        % losses)
        pel_alg_sink_loss_vec = zeros(1,p.Xn-1);
        % looping over the bottom elements
        
        
        b_idx = p.Zn-1;
        for i = 1:p.Xn-1 % looping over bottom elements, starting at the bottom
            pel_alg_sink_loss_vec(i) = p.vA.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2); % [mgC/day]
            b_idx = b_idx + p.Zn-1 -i;
        end
        
        % we divide by the standing stock of the pelagic algae to get the
        % percentage loss per day from sinking
        
        perc_sinking_losses_pel =  sum(pel_alg_sink_loss_vec)/ sum(A.*p.volumes_cyl'); % unit: [1/day]
        
        
        surface_area =  sum(pi*p.W^2);
        
        
        % more response variables
        phyto_vol_avg = (1/p.q)* sum(n_algae_end) ./sum(p.volumes_cyl); % average pelagic algae concentration [mgC/m^3]
        standing_stock = (1/p.q)* ( sum(n_algae_end) + sum(n_benthic_end))./ surface_area; % [mgC/m^2] lake average for the sum of benthic and pelagic algae
        perc_benth = sum(n_benthic_end) /(sum(n_benthic_end) + sum(n_algae_end)); % percentage of total biomass in benthic algae
        pelagic_primary_prod = phytoplankton_primary_prod.*p.volumes_cyl;
        benth_primary_prod = benth_prod.*p.Area_bottom_cyl';
        total_lake_primary_prod = sum( pelagic_primary_prod) + sum(benth_primary_prod);
        
        pelagic_primary_prod_avg_per_area = sum(pelagic_primary_prod)/surface_area;
        benthic_primary_prod_avg_per_area = sum(benth_primary_prod)/surface_area;
        
        
        total_pel_nutrients = sum((n_algae_end + n_dissolved_end + n_detritus_end ).*p.volumes_cyl ); % total ammount of dissolved nutrients in the pelagic;
        
        % light_limited_vol_portion =  sum(p.volumes_cyl(limitation_vec == val*1.0001)) / sum(p.volumes_cyl); % portion of lake that is light limited
        
        
        str  = append( '______________________________________________', '\n');
        str  = append(str, 'Norm of the time derivatives at steady state: ', num2str(x), ' \n ');
        str  = append(str , 'portion of leaked nutrients: ', num2str(abs((p.ntot_0-ntot_end) /p.ntot_0 )), ' \n ' );
        str  = append(str ,'leaked nutrients [mgP]:      ', num2str(abs(p.ntot_0-ntot_end)), ' \n ');
        str  = append(str , '______________________________________________', '\n');
        str  = append(str ,'portion of nutrients in phytoplankton:         ', num2str(perc_algae), ' \n ');
        str  = append(str ,'portion of nutrients in dissolved nutrients: ', num2str(perc_dissolved), ' \n ');
        str  = append(str ,'portion of nutrients in detritus:                   ', num2str(perc_detritus), ' \n ');
        str  = append(str ,'portion of nutrients in sediment:                ', num2str(perc_sediment), ' \n ');
        str  = append(str ,'portion of nutrients in benthic algae:         ', num2str(perc_benthic), '\n');
        str  = append(str , '______________________________________________', '\n');
        str  = append(str ,'portion of pelagic nutrients in pel. algae:         ', num2str(sum((n_algae_end).*p.volumes_cyl)/total_pel_nutrients), '\n');
        str  = append(str ,'portion of pelagic nutrients in dissolved ntr.:   ', num2str(sum((n_dissolved_end).*p.volumes_cyl)/total_pel_nutrients), '\n');
        str  = append(str ,'portion of pelagic nutrients in detritus:             ', num2str(sum((n_detritus_end ).*p.volumes_cyl)/total_pel_nutrients), '\n');
        
        
        dim = [.01 .1 .2 .2];
        textbox = annotation('textbox',dim,'String',sprintf(str),'FitBoxToText','on','EdgeColor','none');
        textbox.FontSize = font_size-4;
        
        
        
        str  = append( '__________________________________________________________', '\n');
        str  = append(str, 'Average pelagic algal biomass [mgC/m^3]:                ', num2str(phyto_vol_avg), '\n');
        str  = append(str, 'Average pelagic dissolved nutrients [mgP/m^3]:         ', num2str(sum(n_dissolved_end)/sum(p.volumes_cyl)), '\n');
        str  = append(str, 'Average pelagic total nutrients [mgP/m^3]:                ', num2str(( sum(n_dissolved_end)+ sum(n_algae_end)+sum(n_detritus_end) )/sum(p.volumes_cyl)), '\n');
        str  = append(str, 'standing stock of pelagic + benth algae [mgC/m^2]:  ', num2str(standing_stock), '\n');
        str  = append(str, 'portion of total biomass in benthic algae :                    ', num2str(perc_benth), '\n');
        
        
        str  = append(str, '__________________________________________________________', '\n');
        str  = append(str, 'lake total primary production [mgC/day]:                       ', num2str(total_lake_primary_prod), '\n');
        str  = append(str, 'benthic contribution to total primary production:            ', num2str(sum(benth_primary_prod)/total_lake_primary_prod), '\n');
        str  = append(str, 'average pelagic primary production [mgC/(day m^2)]:   ', num2str(pelagic_primary_prod_avg_per_area), '\n');
        str  = append(str, 'average benthic primary production [mgC/(day m^2)]:   ', num2str(benthic_primary_prod_avg_per_area), '\n');
        str  = append(str, 'average total primary production [mgC/(day m^2)]:        ', num2str(total_lake_primary_prod/surface_area), '\n');
        
        %str  = append(str, '__________________________________________________________', '\n');
        
        %     str  = append(str, 'portion of lake volume that is light limited: ', num2str(light_limited_vol_portion ), '\n');
        
        dim = [.31 .1 .2 .2];
        textbox = annotation('textbox',dim,'String',sprintf(str),'FitBoxToText','on','EdgeColor','none');
        textbox.FontSize = font_size-4;
        
        
        
        str  = append( '____________________________________________________', '\n');
        str  = append(str, 'portion of light limited pelagic:                       ', num2str( sum(l_lim_vec_pel.*p.volumes_cyl')/sum(p.volumes_cyl)), '\n');
        str  = append(str, 'portion of light limited benthic:                       ', num2str(sum(benthic_limitation.*p.Area_bottom_cyl')/surface_area), '\n');
        
        str  = append(str, '____________________________________________________', '\n');
        str  = append(str, 'pelagic algae sinking losses [1/day]: ', num2str(perc_sinking_losses_pel),'\n');
        dim = [.67 .1 .2 .2];
        textbox = annotation('textbox',dim,'String',sprintf(str),'FitBoxToText','on','EdgeColor','none');
        textbox.FontSize = font_size-4;
        
    end
    
    %% save figure
    
    if(p.model_version == 5)
        fig_name = "2D_results_V5_";
    end
    
    if(p.model_version == 6)
        fig_name = "2D_results_V6_";
    end
    
    
    fig_name = fig_name + "res_" + strrep(num2str(p.Xn),'.','_') + "_kbg_"+ strrep(num2str(p.kbg),'.','_');
    
    fig_name = fig_name + "_max_depth_" + num2str(p.Lmax);
    fig_name = fig_name + "_width_" + num2str(p.W);
    
    % horizontal diffusion coeff (assuming its constant)
    fig_name = fig_name + "_dx_" + strrep(num2str(p.dx(1,1)),'.','_');
    fig_name = fig_name + "_dz_" + strrep(num2str(p.dz(1,1)),'.','_');
    
    if(p.stratified)
        fig_name = fig_name + "_therm_depth_" + strrep(num2str(p.thermocline_depth),'.','_') +  "_therm_thickness_" + strrep(num2str(p.thermocline_thickness),'.','_');
        fig_name = fig_name + "diff_profile_" + num2str(p.diff_above_thermocline) +"," + num2str( p.diff_in_thermocline) + "," + num2str( p.diff_below_thermocline);
    end
    
    
    
    
    if(p.constant_resuspension == 1)
        fig_name = fig_name + "_resus_rate_" + strrep(num2str(p.resus(1)),'.','_') ;
    elseif(p.stepFun_resus)
        fig_name = fig_name +  "_stepFun_resus_upper_th_" + strrep(num2str(p.upper_thresh),'.','_') + "_lower_th_" +strrep(num2str(p.lower_thresh),'.','_') + "_lower_resus_" + strrep(num2str(p.lower_resusp),'.','_')+ "_upper_resus_ " +strrep(num2str(p.upper_resusp),'.','_') ;
        
    else
        fig_name = fig_name +  "_variable_resus_resp_fn_" + num2str(p.response_type);
        
        if(p.manual_in_therm_resus)
            fig_name = fig_name +  "manual_in_therm_resus_" + strrep(num2str(p.manual_therm_resus_val),'.','_');
        end
    end
    
    
    
    
    fig_name = fig_name + '.jpg';
    
    
    exportgraphics(gcf,fig_name);
    %saveas(gcf,fig_name,'epsc');
    
    %% plots of state variable averages for the seasonal models
    
    if(false)
        if(isfield(p,'seasonality'))
            if(p.seasonality)
                figure(figureIndex)
                figureIndex = figureIndex +1;
                
                nutrient_vec = zeros(1,365); % standing stock of dissolved nutrient content over a year [mgP/m^2]
                pelagic_vec = zeros(1,365); % standing stock of pelagic algae content over a year [mgP/m^2]
                benthic_vec = zeros(1,365); % standing stock of bentic algae content over a year [mgP/m^2]
                
                for i=1:365  % looping over the year
                    
                    A = Y_t(t(end-365+i),1:sum(1:p.Xn-1)); % phytoplankton
                    Rd = Y_t(t(end-365+i),sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
                    B = Y_t(t(end-365+i), 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae
                    
                    
                    pelagic_vec(i) = sum(A'.*p.q.*p.volumes_cyl)/sum(p.Area_bottom_cyl);
                    nutrient_vec(i) = sum(Rd'.*p.volumes_cyl)/sum(p.Area_bottom_cyl);
                    benthic_vec(i) = sum(B.*p.q.*p.Area_bottom_cyl)/sum(p.Area_bottom_cyl);
                    
                end
                
                plot([1:365],nutrient_vec,'color','black','linewidth',1.5);
                hold on
                plot([1:365],pelagic_vec,'color','green','linewidth',1.5);
                hold on
                plot([1:365],benthic_vec,'color','blue','linewidth',1.5);
                title("Standing stock [mgP/m^2]. kbg = " +num2str(p.kbg)+ ", resus type = " + num2str(p.response_type));
                legend('dissolved nutrients','pelagic algae','benthic algae')
                xlabel('day of the year');
                ylabel('standing stock [mgP/m^2]');
                
                
                % save figure
                
                fig_name = "standing_stock_V6_" + "res_" + strrep(num2str(p.Xn),'.','_') + "_kbg_"+ strrep(num2str(p.kbg),'.','_');
                
                if(p.stratified)
                    if(~p.seasonality_thermoC)
                        fig_name = fig_name + "_therm_depth_" + strrep(num2str(p.thermocline_depth),'.','_') +  "_therm_thickness_" + strrep(num2str(p.thermocline_thickness),'.','_');
                    end
                    fig_name = fig_name + "_diff_profile_" + num2str(p.diff_above_thermocline) +"," + num2str( p.diff_in_thermocline) + "," + num2str( p.diff_below_thermocline);
                end
                
                if(p.constant_resuspension == 1)
                    fig_name = fig_name + "_resus_rate_" + strrep(num2str(p.resus(1)),'.','_') ;
                elseif(p.stepFun_resus)
                    fig_name = fig_name +  "_stepFun_resus_upper_th_" + strrep(num2str(p.upper_thresh),'.','_') + "_lower_th_" +strrep(num2str(p.lower_thresh),'.','_') + "_lower_resus_" + strrep(num2str(p.lower_resusp),'.','_')+ "_upper_resus_ " +strrep(num2str(p.upper_resusp),'.','_') ;
                    
                else
                    fig_name = fig_name +  "_variable_resus_resp_fn_" + num2str(p.response_type);
                    
                    if(p.manual_in_therm_resus)
                        fig_name = fig_name +  "manual_in_therm_resus_" + strrep(num2str(p.manual_therm_resus_val),'.','_');
                    end
                end
                
                
                saveas(gcf,fig_name,'epsc');
                
            end
        end
    end
    
    %% cross radial profiles
    %     figure(figureIndex)
    %     figureIndex = figureIndex +1;
    %
    %     % calculating the light/nutrient limitation for the benthic algae
    %     b_idx = 1;
    %     b_idx_next = p.Zn-1;
    %     depth_integrated_pel_algae = zeros(1,p.Xn-1); % depth integrated algal nutrient content [mgP/m^2]
    %     benthic_limitation = zeros(p.Xn-1,1);
    %
    %     for i=1:p.Xn-1
    %
    %
    %         depth_integrated_pel_algae(i) = (sum(A(b_idx:b_idx_next).*p.volumes_cyl(b_idx:b_idx_next)'))/p.Area_bottom_cyl(i) ;
    %
    %         size(b_idx:b_idx_next)
    %         b_idx = b_idx_next+1;
    %         b_idx_next = b_idx-1 + p.Zn-1 -i;
    %     end
    
end


test = 1;

end