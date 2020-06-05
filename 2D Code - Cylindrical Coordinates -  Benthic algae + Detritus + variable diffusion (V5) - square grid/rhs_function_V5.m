function [output] = rhs_function_V5(t,Y,p)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface layer and so forth.

% function is used as input to ode15s.

% Version 5 - detritus + variable diffusion coefficients. cut rectangular
% grid


%% reformatting input
% separating the state variables from y

A = Y(1:sum(1:p.Xn-1)); % phytoplankton
Rd = Y(sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
D = Y(2*sum(1:p.Xn-1) +1 : 3*sum(1:p.Xn-1)); % ditritus
Rs = Y(3*sum(1:p.Xn-1)+1 : 3*sum(1:p.Xn-1)+ p.Xn-1); % sedimented nutrients
B = Y( 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae

%% Calculation of the light intensity at each grid element center
%A_temp = A';
%D_temp = D';
%Z_vol_temp = p.Z_vol';

I = zeros(sum(1:p.Xn-1),1);

index = 1;
temp = 1;

for j = 1:p.Xn-1 % pooling over each column (j = column nr)
    int = 0;
    for i = 1:p.Zn - temp % traversing down one column (i = row nr)
        
        dz_middle = p.Z_vol(i,j) - p.Z(i,j); % distance from bottom of grid element (i-1,j) to the center of grid element (i,j)
        
        if(i==1) % width of grid element (i,j) in the z-direction
            dz = p.Z(i+1,j);
        else
            dz = p.Z(i+1,j) - p.Z(i,j);
        end
        
        int_middle = int + (p.kD.*1/p.q.*D(index) +p.kA*A(index)+ p.kbg)*dz_middle; % this is the integral down to the center of grid element (i,j).
        int = int + (p.kD.*1/p.q.*D(index) +p.kA*A(index)+ p.kbg)*dz; % here we integrate to the bottom of grid element (i,j)
        I(index) =  p.I0 .* exp( - int_middle);
        index = index +1;
    end
    temp = temp +1;
end
p.I = I;

%% Lake bulk dynamics
dAdt = zeros(sum(1:p.Xn-1),1);
dRdt = zeros(sum(1:p.Xn-1),1);
dDdt = zeros(sum(1:p.Xn-1),1);
dRsdt = zeros(1,p.Xn-1);
dBdt = zeros(1,p.Xn-1);

% Growth, death, and respiration losses of the algae
G = p.Gmax .* min(Rd./(Rd+p.M), I./(I+p.H));
dAdt = dAdt + (G- p.Ad - p.lbg_A).*A;

% Dead algae wind up as detritus in the water, and detritus remineralize
dDdt = dDdt +  p.q*p.Ad.*A - p.Dbg*D;

% Alage consume and respire nutrients, nutrients from detritus remineralize
dRdt =  dRdt +  p.q.*(p.lbg_A-G).*A + p.Dbg*D;

%% Diffusion & sinking

% reformatting bulk variables into matrixes
A_matrix  = NaN.*ones(p.Zn-1, p.Xn-1);
Rd_matrix = NaN.*ones(p.Zn-1, p.Xn-1);
D_matrix  = NaN.*ones(p.Zn-1, p.Xn-1);

temp = 1;
index = 1;
for j = 1:p.Xn-1 % looping over each column
    for i = 1:p.Zn - temp % traversing down one column (i = row nr)
        %I_matrix(i,j) = I(index);
        A_matrix(i,j)  = A(index);
        Rd_matrix(i,j) = Rd(index);
        D_matrix(i,j)  = D(index);
        index = index +1;
    end
    temp = temp +1;
end

dAdt_matrix = zeros(length(A),1);
dRdt_matrix = zeros(length(A),1);
dDdt_matrix = zeros(length(A),1);

index =1;

for j = 1:p.Xn-1
    for i = 1:p.Zn-j
        
        % diffusion & sinking
        
        if(i>1) % top
            
            % diffusion
            temp_A =  0.5*(p.X(i,j)^2 - p.X(i,j+1)^2)*0.5*(p.dz(i,j) + p.dz(i-1,j))*(A_matrix(i,j) - A_matrix(i-1,j))/(p.Z_vol(i,j)-p.Z_vol(i-1,j));
            temp_Rd = 0.5*(p.X(i,j)^2 - p.X(i,j+1)^2)*0.5*(p.dz(i,j) + p.dz(i-1,j))*(Rd_matrix(i,j) - Rd_matrix(i-1,j))/(p.Z_vol(i,j)-p.Z_vol(i-1,j));
            temp_D =  0.5*(p.X(i,j)^2 - p.X(i,j+1)^2)*0.5*(p.dz(i,j) + p.dz(i-1,j))*(D_matrix(i,j) - D_matrix(i-1,j))/(p.Z_vol(i,j)-p.Z_vol(i-1,j));
                        
            if(~isnan(temp_A))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A;
            end            
            if(~isnan(temp_Rd))
                dRdt_matrix(index) = dRdt_matrix(index) + temp_Rd;
            end            
            if(~isnan(temp_D))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D;
            end
            
            % sinking
            temp_A_sink = 0.5*(p.X(i,j+1)^2 - p.X(i,j)^2)*p.vA*0.5*(A_matrix(i,j) + A_matrix(i-1,j));
            temp_D_sink = 0.5*(p.X(i,j+1)^2 - p.X(i,j)^2)*p.vD*0.5*(D_matrix(i,j) + D_matrix(i-1,j));
            
             if(~isnan(temp_A_sink))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A_sink;
            end            
            if(~isnan(temp_D_sink))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D_sink;
            end            
            
            
        end
        
        if(i<p.Zn-1) % bottom integral
            
            temp_A  = 0.5*(p.X(i+1,j+1)^2-p.X(i+1,j)^2)*0.5*(p.dz(i,j)+p.dz(i+1,j))*(A_matrix(i+1,j)- A_matrix(i,j))/(p.Z_vol(i+1,j)-p.Z_vol(i,j));
            temp_Rd = 0.5*(p.X(i+1,j+1)^2-p.X(i+1,j)^2)*0.5*(p.dz(i,j)+p.dz(i+1,j))*(Rd_matrix(i+1,j)- Rd_matrix(i,j))/(p.Z_vol(i+1,j)-p.Z_vol(i,j));
            temp_D  = 0.5*(p.X(i+1,j+1)^2-p.X(i+1,j)^2)*0.5*(p.dz(i,j)+p.dz(i+1,j))*(D_matrix(i+1,j)- D_matrix(i,j))/(p.Z_vol(i+1,j)-p.Z_vol(i,j));
            
            if(~isnan(temp_A))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A;
            end            
            if(~isnan(temp_Rd))
                dRdt_matrix(index) = dRdt_matrix(index) + temp_Rd;
            end            
            if(~isnan(temp_D))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D;
            end
            
            % sinking
            temp_A_sink = -0.5*(p.X(i+1,j+1)^2-p.X(i+1,j)^2)*p.vA*0.5*(A_matrix(i,j) + A_matrix(i+1,j));
            temp_D_sink = -0.5*(p.X(i+1,j+1)^2-p.X(i+1,j)^2)*p.vD*0.5*(D_matrix(i,j) + D_matrix(i+1,j));
            
             if(~isnan(temp_A_sink))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A_sink;
            end            
            if(~isnan(temp_D_sink))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D_sink;
            end               
        end
        
        if(j>1) % left integral
           temp_A =  -p.X(i,j)*(p.Z(i+1,j)-p.Z(i,j))*0.5*(p.dx(i,j)+p.dx(i,j-1))*(A_matrix(i,j)- A_matrix(i,j-1))/(p.X_vol(i,j)-p.X_vol(i,j-1));
           temp_Rd = -p.X(i,j)*(p.Z(i+1,j)-p.Z(i,j))*0.5*(p.dx(i,j)+p.dx(i,j-1))*(Rd_matrix(i,j)- Rd_matrix(i,j-1))/(p.X_vol(i,j)-p.X_vol(i,j-1));
           temp_D =  -p.X(i,j)*(p.Z(i+1,j)-p.Z(i,j))*0.5*(p.dx(i,j)+p.dx(i,j-1))*(D_matrix(i,j)- D_matrix(i,j-1))/(p.X_vol(i,j)-p.X_vol(i,j-1));
            
            if(~isnan(temp_A))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A;
            end            
            if(~isnan(temp_Rd))
                dRdt_matrix(index) = dRdt_matrix(index) + temp_Rd;
            end            
            if(~isnan(temp_D))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D;
            end
        end
        
        if(j<p.Xn-1) % right integral
            
            temp_A  = p.X(i,j+1)*(p.Z(i+1,j+1)-p.Z(i,j+1))*0.5*(p.dx(i,j)+p.dx(i,j+1))*(A_matrix(i,j+1)- A_matrix(i,j))/(p.X_vol(i,j+1)-p.X_vol(i,j));
            temp_Rd = p.X(i,j+1)*(p.Z(i+1,j+1)-p.Z(i,j+1))*0.5*(p.dx(i,j)+p.dx(i,j+1))*(Rd_matrix(i,j+1)- Rd_matrix(i,j))/(p.X_vol(i,j+1)-p.X_vol(i,j));
            temp_D  = p.X(i,j+1)*(p.Z(i+1,j+1)-p.Z(i,j+1))*0.5*(p.dx(i,j)+p.dx(i,j+1))*(D_matrix(i,j+1)- D_matrix(i,j))/(p.X_vol(i,j+1)-p.X_vol(i,j));
            
            if(~isnan(temp_A))
                dAdt_matrix(index) = dAdt_matrix(index) + temp_A;
            end            
            if(~isnan(temp_Rd))
                dRdt_matrix(index) = dRdt_matrix(index) + temp_Rd;
            end            
            if(~isnan(temp_D))
                dDdt_matrix(index) = dDdt_matrix(index) + temp_D;
            end
        end
        
        
        
        index = index +1;
        
        %sinking
        
    end
end

% Division by element volume (from l.h.s. when integrating over a control
% volume). The factor 2pi comes from the radial integration, and since we
% assume symmetry the integrand is independent of r and we can simply
% integrate the radial integral to 2pi.

dAdt_matrix = dAdt_matrix.*2.*pi./p.volumes_cyl; 
dDdt_matrix = dDdt_matrix.*2.*pi./p.volumes_cyl;
dRdt_matrix = dRdt_matrix.*2.*pi./p.volumes_cyl;


dAdt = dAdt + dAdt_matrix;
dRdt = dRdt + dRdt_matrix;
dDdt = dDdt + dDdt_matrix;

%% Bottom dynamics
if(true)
    % looping over the bottom elements
    b_idx = p.Zn-1; % bottom index
    for i = 1:p.Xn-1 % looping over bottom elements, starting at the bottom
        
        % Bethic algae net growth
        nutrient_limited_growth =  p.Gmax_benth.*B(i)'.*(Rd(b_idx)./(Rd(b_idx)+p.M_benth));
        light_limited_growth = p.Gmax_benth./p.kB.*log((p.H_benth + I(b_idx))./(p.H_benth + I(b_idx).*exp(-p.kB.*B(i))));
        
        dBdt(i) = dBdt(i) + min(nutrient_limited_growth, light_limited_growth) - p.lbg_benth*B(i);
        
        % Nutrients in the sediment remineralize into the lake.
        dRsdt(i) =  dRsdt(i) - p.r.*Rs(i);
        
        % Benthic Algae respire/die and a portion p.benth_recycling is bound in
        % particulate form and is introduced into the sediment. 
        dRsdt(i) =  dRsdt(i) + (1-p.benth_recycling).*p.q_benth.*p.lbg_benth*B(i);
        
        % Algae sinks into the sediment, dies and instantly becomes sedimented nutrients.        
        dAdt(b_idx) = dAdt(b_idx) - p.vA.*p.death_rate.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 -  p.X(end,i)^2)./p.volumes_cyl(b_idx);
           dRsdt(i) = dRsdt(i) + p.vA.*p.death_rate.*p.q.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.Area_bottom_cyl(i);
        
        % Detritus sinks into the sediment, dies and instantly becomes sedimented  nutrients.
       dDdt(b_idx) = dDdt(b_idx) - p.vD.*D(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.volumes_cyl(b_idx);
          dRsdt(i) = dRsdt(i) + p.vD.*D(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.Area_bottom_cyl(i);
        
        
        % Detritus is resuspended into the water
        dRsdt(i) =  dRsdt(i) - p.resus.*Rs(i);
        dDdt(b_idx) = dDdt(b_idx) +  p.resus.*Rs(i).*p.Area_bottom_cyl(i)/p.volumes_cyl(b_idx);
        
        
        %%%%% net flux of dissolved nutrients into the lake from the sediment and consumption by the benthic algae   %%%%%%%%%%%%%%%%%%%%%%%
        
        nutrient_limited_flux =  p.q_benth.*B(i)'.*p.Gmax_benth.*(Rd(b_idx)./(Rd(b_idx)+p.M_benth)); % [mgP/ (m^2 day)]
        light_limited_flux = p.q_benth.*p.Gmax./p.kB*log((p.H_benth + I(b_idx))./(p.H_benth +I(b_idx).*exp(-p.kB.*B(i)'))); % [mgP/ (m^2 day)]
        
        benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,
        net_flux = p.r.*Rs(i) +  p.benth_recycling.*p.q_benth.*p.lbg_benth*B(i)' - benth_P_consumption;
        
        dRdt(b_idx) =   dRdt(b_idx) + net_flux.*p.Area_bottom_cyl(i)./p.volumes_cyl(b_idx); % the net flux here has the units [mg P/ (m^2 day)], and is multiplied by the area of the bttom segment to yield
        % a change [mg P/ day]. Division by the element volume yields the sought change in concentration.
        
        b_idx = b_idx + p.Zn-1 -i;
    end
end

%% Output

output = [dAdt(:); dRdt(:); dDdt(:); dRsdt(:); dBdt(:)]; % The output is a column vector with the derivatives of the state variables at each grid point.


end
