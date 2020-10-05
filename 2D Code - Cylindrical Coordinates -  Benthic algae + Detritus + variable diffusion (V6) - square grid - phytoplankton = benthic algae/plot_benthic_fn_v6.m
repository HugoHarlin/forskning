
function [] = plot_benthic_fn_v6(t,Y_t,p)
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


%% Benthic algae
figure(1)
font_size = 16;
max_survival_depth = -1./p.kbg .* log(p.H_benth./( p.I0.*(p.Gmax_benth./p.lbg_benth -1) )); % maximum survival depth of benthic algae assuming no phytoplankton
  max_survival_depth_x = p.W*((max_survival_depth./p.Lmax -1)*(p.Lmax/(p.Lmin-p.Lmax)))^(1/p.alpha); % corressponding distance from the lake center
     
sed_points = zeros(1,p.Xn-1); % creating points for sediment plot
sed_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
query_points = linspace(0,p.W,200);
%sed_interp = interp1(sed_points,Rs, query_points,'makima');
sed_interp = interp1(sed_points,Rs, query_points,'pchip');


% creating points for benthic algae plot
benth_points = zeros(1,p.Xn-1);
benth_points(:) = (p.X(1,2:end) + p.X(1,1:end-1))/2;
new_res = 300; % resolution of interpolated plot
query_points = linspace(0,p.W,new_res);
%benth_interp = interp1(sed_points,B, query_points,'makima');
benth_interp = interp1(sed_points,B.*p.q, query_points,'pchip','extrap'); % interpolated net growtht values, unit [mgP/(m^2 day)]


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
  %  if((nutrient_limited_growth(i) < 0)) nutrient_limited_growth(i) = 0; end
  %  if((light_limited_growth(i) < 0)) light_limited_growth(i) = 0; end
    
    
    
    %   benth_prod_interp(i) = min(nutrient_limited_growth(i),light_limited_growth(i))/Benth_interp(i);
    
    
    
    if(Benth_interp(i) > 1e-8) % threshold value, below which the benthic algae is considered extinct and we don't compute the growth limitation (which requires algae present. if not done the specific growth plot looks very erratic where there is miniscule/no ammounts of benthic algae.
        benth_prod_interp(i) = min(nutrient_limited_growth(i),light_limited_growth(i))/Benth_interp(i);
    else
       % benth_prod_interp(i) = 0; % the calculation of the benthic growth limitation requires there to be benthic algae present since the benthic algae shadows itself and this is taken into account.
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
ylabel("concentration [mgP/m^2]");
xlabel("distance from lake center [m]");
ylim([0 inf]);


% we add a second plot of the specific production on top of the
% benthic nutrient content, adding a second y-axis to the right
hold on

% right axis
yyaxis right
test = plot(query_points, benth_prod_interp,'color','black','lineWidth',2);
set(test,'Linestyle','-'); 
ylabel("benthic specific production");
ylabel("specific production");

set(gca,'fontSize',font_size);
pbaspect([1.2 1 1])

   %% save figure
    
    if(p.model_version == 5)
        fig_name = "2D_benthic_V5_";
    end
    
    if(p.model_version == 6)
        fig_name = "2D_benthic_V6_";
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


end
