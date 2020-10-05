function [output] = rhs_function_V6(t,Y,p)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface layer and so forth.

% function is used as input to ode15s.

% Version 6 - detritus + variable diffusion coefficients. cut rectangular
% grid

%t
%% reformatting input
% separating the state variables from y

A = Y(1:sum(1:p.Xn-1)); % phytoplankton
Rd = Y(sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
D = Y(2*sum(1:p.Xn-1) +1 : 3*sum(1:p.Xn-1)); % ditritus
Rs = Y(3*sum(1:p.Xn-1)+1 : 3*sum(1:p.Xn-1)+ p.Xn-1); % sedimented nutrients
B = Y( 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae

%% Light intensity and seasonality
%A_temp = A';
%D_temp = D';
%Z_vol_temp = p.Z_vol';

I = zeros(sum(1:p.Xn-1),1);

index = 1;
temp = 1;


% if p.seasonality = 1 then seasonality is activated.
% I added additional flags for the seasonally varying parameters so
% these can be tracked after a simulation is done (i.e.
% p.seasonality_light need to be true in addition for p.seasonality for a
% seasonally varying light intesity).

if(isfield(p,'seasonality'))
    if(p.seasonality)
        if(p.seasonality_light)
            % seasonal light is modelled with a sin-function that varies from
            % p.minLight to p.maxlight
            p.I0 = (p.maxLight-p.minLight)/2*sin(2*pi*(90+mod(t(end),365))/365) + (p.minLight + p.maxLight)/2;
        end
        
        
        
        if(p.stratified)
            if(isfield(p,'seasonality_thermoC'))
                if(p.seasonality_thermoC) % varying depth of the thermocline
                    
                    
                    if(false) % sinusoidal anuual variation of thermocline depth
                        p.thermocline_depth = (p.maxTherm -p.minTherm)/2*sin(2*pi*(-90+mod(t(end),365))/365) + (p.minTherm + p.maxTherm )/2;
                    end
                    
                    
                    if(true) % manually constructed thermocline seasonality.
                        dx = ones(p.Zn-1,p.Xn-1); % Radial Turbulent-diffusion coefficient [m^2 day^-1]
                        dz = ones(p.Zn-1,p.Xn-1); % Vertical Turbulent-diffusion coefficient [m^2 day^-1]
                        
                        
                        if(mod(t(end),365)< 61 && mod(t(end),365) > 29 ) % days 30-60 are considered winter years, the thermocline is not present and the diffusion is the same everywhere in the lake.
                            dx = p.diff_above_thermocline.*dx;
                            dz = p.diff_above_thermocline.*dz;
                            
                            p.dx = dx;
                            p.dz = dz;
                        else % we are not in the winter months and the lake is stratified.
                            
                            if(mod(t(end),365)> 240 ) % We are not in the summer period and the thermocline depth is linearly increasing
                                p.thermocline_depth = p.minTherm + (p.maxTherm - p.minTherm)*(mod(t(end),365)-240)/(155);
                            elseif(mod(t(end),365)<= 29)
                                % start of the year before well mixed winter
                                % period.
                                p.thermocline_depth = (p.minTherm + (p.maxTherm - p.minTherm)*(364-240)/(155)) +(p.maxTherm-(p.minTherm + (p.maxTherm - p.minTherm)*(364-240)/(155)))* (mod(t(end),365))/(29) ;
                            else
                                % summer months, thermocline is at its shallowest.
                                p.thermocline_depth = p.minTherm;
                                
                            end
                            
                            
                            %%%%% manual setting of the diffusion coefficients %%%%
                            for x = 1:p.Xn-1
                                for z = 1:p.Zn-1
                                    
                                    if(p.Z_vol(z,x) < p.thermocline_depth)
                                        p.dx(z,x) =  p.diff_above_thermocline*(0.8 + 0.2*( pp.thermocline_depth- abs(pp.Z_vol(z,x)))/pp.thermocline_depth); % *(p.Zn-z)/p.Zn;
                                        p.dz(z,x) =  p.diff_above_thermocline*(0.8 + 0.2*( pp.thermocline_depth- abs(pp.Z_vol(z,x)))/pp.thermocline_depth); % *(p.Zn-z)/p.Zn;
                                    elseif(p.thermocline_depth < p.Z_vol(z,x) &&  p.Z(z,x) < p.thermocline_depth + p.thermocline_thickness)
                                        p.dx(z,x) = p.diff_in_thermocline;
                                        p.dz(z,x) = p.diff_in_thermocline;
                                    else
                                        p.dx(z,x) =  p.diff_below_thermocline*(0.8 + 0.2*( pp.Lmax- abs(pp.Z_vol(z,x)))/pp.Lmax); % *(p.Zn-0.5*z)/(p.Zn);
                                        p.dz(z,x) =  p.diff_below_thermocline*(0.8 + 0.2*( pp.Lmax- abs(pp.Z_vol(z,x)))/pp.Lmax); % *(p.Zn-0.5*z)/(p.Zn);
                                    end
                                end
                            end
                        end
                        
                    end
                    
                    % removing values outside the grid
                    
                    jj = p.Xn-1;
                    for ii=2:p.Zn-1
                        p.dx(ii,jj:end) = NaN;
                        p.dz(ii,jj:end) = NaN;
                        jj = jj -1;
                    end
                    
                end
            end
            
            
            
            
            % if the thermocline is moved and the resuspension depends on the
            % diffusion coefficients, we need to calculate the resuspension for the
            % new conditions
            
            if(~p.constant_resuspension)
                %%%% Resuspension rate that varies with depth %%%%
                p.resus = ones(1,p.Xn-1);
                
                
                % Type 1 functional response (linear)
                if(p.response_type == 1)
                    max_resus = 0.2;
                    min_resus = 0;
                    b_idx = p.Zn-1; % bottom index
                    for i=1:p.Xn-1
                        %   diff_normal = (p.dx(i,b_idx) + p.dz(i,b_idx)* p.W/p.Lmax)/sqrt(1 + p.W^2/p.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                        diff_normal = p.dz(b_idx,i);
                        p.resus(i) = min_resus +  ((max_resus - min_resus)/100)*diff_normal;
                        b_idx = p.Zn-1 -i;
                    end
                end
                
                % Type 2 functional response
                if(p.response_type == 2)
                    % the parameters are tweaked manually to get a
                    % functional response that had a "good" s-shape to
                    % it, assuming a maimum resuspension rate of 20%/day.
                    
                    c =1; % 4.2;
                    b = 0.3; %2.0e-05;
                    
                    b_idx = p.Zn-1; % bottom index
                    for i=1:p.Xn-1
                        %   diff_normal = (p.dx(i,b_idx) + p.dz(i,b_idx)* p.W/p.Lmax)/sqrt(1 + p.W^2/p.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                        diff_normal = p.dz(b_idx,i);
                        p.resus(i) = (0.2*b.*diff_normal.^c)./(1 + 1*b.*diff_normal.^c);
                        b_idx = p.Zn-1 -i;
                    end
                    test = 1;
                end
                
                
                % Type 3 functional response
                if(p.response_type == 3)
                    % the parameters are tweaked manually to get a
                    % functional response that had a "good" s-shape to
                    % it, assuming a maimum resuspension rate of 20%/day.
                    
                    c =3; % 4.2;
                    b = 0.02; %2.0e-05;
                    
                    b_idx = p.Zn-1; % bottom index
                    for i=1:p.Xn-1
                        %   diff_normal = (p.dx(i,b_idx) + p.dz(i,b_idx)* p.W/p.Lmax)/sqrt(1 + p.W^2/p.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                        diff_normal = p.dz(b_idx,i);
                        p.resus(i) = (0.2*b.*diff_normal.^c)./(100 + 1*b.*diff_normal.^c);
                        b_idx = p.Zn-1 -i;
                    end
                end
                
                % Type 4 functional response
                if(p.response_type == 4)
                    
                    % Dependent on bottom shear stress tau, which is calculated by
                    % assuming that close to the lake floor there is only movement
                    % parallell to the bottom, and the change in concentration due
                    % to diffusion is interpreted as a fictive speed
                    % v_fic which is parallel to the lake floor. in mathematical terms:
                    
                    % Dr*d^2A/dr^2 = v_fic * dA/dr, (1)
                    
                    % where rhat is parallell to the
                    % bottom.
                    
                    % Solving for v_fic we get an expression for the speed of the
                    % Algae (and the surrounding water causing the resuspension)
                    % which we can use to calculate the shear stress on the
                    % bottom.
                    
                    % The question is if equation (1) makes sense. The formula we
                    % use to calculate the shear stress assumes a no-slip
                    % condition on the bottom interface, are we breaking this
                    % assumption with equation (1)?
                    
                    % source: [Jin 2004 – "Case study: modeling of sediment
                    % transport and wind-wave impact in lake Okechobee"]
                    
                    tau_thresh = 0.01; % minimum threshold shear stress for resuspension to occur
                    tau = zeros(1,p.Xn-1); % bottom shear stress [N/m^2];
                    b_idx = p.Zn-1; % bottom index
                    
                    % Distance between two adjacent nodes along the lake
                    % bottom.
                    node_dist =  sqrt( (p.X_vol(p.Zn-i-1,i+1) - p.X_vol(p.Zn-i,i))^2 +  (p.Z_vol(p.Zn-i-1,i+1) - p.Z_vol(p.Zn-i,i))^2 );
                    
                    
                    % calculating the resuspension rate at each bottom grid
                    % cell.
                    for i=1:p.Xn-1
                        
                        
                        % first and second derivatives of the phytoplankton
                        % concentration parallell to the bottom.
                        
                        % current node index b_idx{0} = b_idx
                        
                        % node index of previous node (to the left, further down into the
                        % lake towards the middle)
                        % b_idx_{-1} = b_idx-p.Zn+i;
                        
                        
                        % node index of the next node (to the right, furhter
                        % up towards the lake edge.)
                        % b_idx_{+1} = b_idx + p.Zn-1 -i;
                        
                        % and so forth:
                        %  b_idx_{+2} = b_idx + p.Zn-1 -i + p.Zn-1 -(i+1) =  b_idx + 2*(p.Zn-1) -2*i +1
                        %  b_idx_{-2} = b_idx - p.Zn-1 +i - p.Zn-1 + (i-1)  =  b_idx - 2*(p.Zn-1) + 2*i -1
                        
                        
                        if(i==1)
                            dAdr   = (A(b_idx+p.Zn-1-i)-A(b_idx))/node_dist;
                            d2Adr2 = (A(b_idx) -2*A(b_idx + p.Zn-1 -i) + A(b_idx + 2*(p.Zn-1) -2*i +1))/node_dist^2;
                        elseif(i == p.Xn-1)
                            dAdr   = (A(b_idx)-A(b_idx-p.Zn+i))/node_dist;
                            d2Adr2 =  (A(b_idx) -2*A(b_idx-p.Zn+i) +  A(b_idx - 2*(p.Zn-1) + 2*i -1))/node_dist^2;
                        else
                            dAdr   =  0.5*(A(b_idx+p.Zn-1-i)-A(b_idx-p.Zn+i))/node_dist; % central difference 1 order
                            d2Adr2 =  (A(b_idx -(p.Zn-i) )  -2*A(b_idx) + A(b_idx+p.Zn-1-i)) / node_dist^2; % central difference 2nd order
                        end
                        
                        
                        v = p.dx*d2Adr2/dAdr; % fictive water flow speed along the bottom
                        
                        rhoe = 1000;
                        tau = v^2*rhoe;
                        
                        tau_crit = 0.216; % Critical shear stress for erosion [N/m^2]
                        resus_rate = 0.06*1000; % resuspension rate [mgP/m^2]
                        p.resus(i) = 1;
                        b_idx = b_idx + p.Zn-1 -i;
                    end
                end
                
                
                if(p.manual_in_therm_resus) % sets the resuspension coefficient in the thermocline manually
                    b_idx = p.Zn-1; % bottom index
                    for i=1:p.Xn-1
                        %   diff_normal = (p.dx(i,b_idx) + p.dz(i,b_idx)* p.W/p.Lmax)/sqrt(1 + p.W^2/p.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                        if(p.dz(b_idx,i) == p.diff_in_thermocline)
                            p.resus(i) = p.manual_therm_resus_val;
                        end
                        b_idx = p.Zn-1 -i;
                    end
                    
                end
                
                
                % stepfunction resuspension
                if(p.stepFun_resus)
                    
                    p.upper_resusp = 0.05; % resuspension rate above upper threshold [1/day]
                    p.lower_resusp = 0.0001; % resuspension rate below lower threshold [1/day]
                    
                    p.upper_thresh = 10; % upper threshold depth [m]
                    p.lower_thresh = 12; % lower threshold depth [m]
                    
                    
                    for i=1:p.Xn-1
                        if(p.Z_vol(p.Xn-i) < p.upper_thresh) % above upper threshold.
                            p.resus(i) = p.upper_resusp;
                            
                            
                        elseif(p.Z_vol(p.Xn-i) > p.lower_thresh) % below lower threshold.
                            p.resus(i) = p.lower_resusp;
                            
                        else % we are between the thresholds
                            
                            % diff_normal = (p.dx(i,b_idx) + p.dz(i,b_idx)* p.W/p.Lmax)/sqrt(1 + p.W^2/p.Lmax^2); % diffusion coefficient in the normal direction with respect to the (linear) bottom.
                            % diff_normal = p.dz(b_idx,i);
                            p.resus(i) =   p.lower_resusp +    (p.uper_resusp - p.lower_resusp)*(1- ( p.Z_vol(p.Xn-i)- p.uper_thresh )/(p.lower_thresh - p.upper_thresh));
                            
                            
                        end
                    end
                end
                
                p.resus_benth = p.resus;
            end
            
        end
        
    end
end


% the resuspension code below is duplicated from the seasonality bit above.
% I wanted to test the type 4 response without seasonality, thats why i
% copied it here.
if(~p.constant_resuspension)
    
    % Type 4 functional response
    if(p.response_type == 4)
        
        % Dependent on bottom shear stress tau, which is calculated by
        % assuming that close to the lake floor there is only movement
        % parallell to the bottom, and the change in concentration due
        % to diffusion is interpreted as a fictive speed
        % v_fic which is parallel to the lake floor. in mathematical terms:
        
        % Dr*d^2A/dr^2 = v_fic * dA/dr, (1)
        
        % where rhat is parallell to the
        % bottom.
        
        % Solving for v_fic we get an expression for the speed of the
        % Algae (and the surrounding water causing the resuspension)
        % which we can use to calculate the shear stress on the
        % bottom.
        
        % The question is if equation (1) makes sense. The formula we
        % use to calculate the shear stress assumes a no-slip
        % condition on the bottom interface, are we breaking this
        % assumption with equation (1)?
        
        % source: [Jin 2004 – "Case study: modeling of sediment
        % transport and wind-wave impact in lake Okechobee"]
        
        tau_thresh = 0.01; % minimum threshold shear stress for resuspension to occur
        tau_vec = zeros(1,p.Xn-1); % bottom shear stress [N/m^2];
        b_idx = p.Zn-1; % bottom index
        
        % Distance between two adjacent nodes along the lake
        % bottom.
        i=1;
        node_dist =  sqrt( (p.X_vol(p.Zn-i-1,i+1) - p.X_vol(p.Zn-i,i))^2 +  (p.Z_vol(p.Zn-i-1,i+1) - p.Z_vol(p.Zn-i,i))^2 );
        
        
        % calculating the resuspension rate at each bottom grid
        % cell.
        for i=1:p.Xn-1
            
            
            % first and second derivatives of the phytoplankton
            % concentration parallell to the bottom.
            
            % current node index b_idx{0} = b_idx
            
            % node index of previous node (to the left, further down into the
            % lake towards the middle)
            % b_idx_{-1} = b_idx-p.Zn+i;
            
            
            % node index of the next node (to the right, furhter
            % up towards the lake edge.)
            % b_idx_{+1} = b_idx + p.Zn-1 -i;
            
            % and so forth:
            %  b_idx_{+2} = b_idx + p.Zn-1 -i + p.Zn-1 -(i+1) =  b_idx + 2*(p.Zn-1) -2*i +1
            %  b_idx_{-2} = b_idx - p.Zn-1 +i - p.Zn-1 + (i-1)  =  b_idx - 2*(p.Zn-1) + 2*i -1
            
            
            if(i==1)
                dAdr   = (A(b_idx+p.Zn-1-i)-A(b_idx))/node_dist;
                d2Adr2 = (A(b_idx) -2*A(b_idx + p.Zn-1 -i) + A(b_idx + 2*(p.Zn-1) -2*i +1))/node_dist^2;
            elseif(i == p.Xn-1)
                dAdr   = (A(b_idx)-A(b_idx-p.Zn+i))/node_dist;
                d2Adr2 =  (A(b_idx) -2*A(b_idx-p.Zn+i) +  A(b_idx - 2*(p.Zn-1) + 2*i -1))/node_dist^2;
            else
                dAdr   =  0.5*(A(b_idx+p.Zn-1-i)-A(b_idx-p.Zn+i))/node_dist; % central difference 1 order
                d2Adr2 =  (A(b_idx -(p.Zn-i) )  -2*A(b_idx) + A(b_idx+p.Zn-1-i)) / node_dist^2; % central difference 2nd order
            end
            
            
            v = p.dx*d2Adr2/dAdr; % fictive water flow speed along the bottom
            
            rhoe = 1000;
            tau = v^2*rhoe;
            
            tau_crit = 0.216; % Critical shear stress for erosion [N/m^2]
            resus_rate = 0.06*1000; % resuspension rate [mgP/m^2] (isn't this is what I'm seeking, what is the difference between the resus rate and the resus flux?)
            alpha = 1; % exponent in resuspension formula
            
            
            if(tau >= tau_crit) % if the bed stress is higher than the critical bed stress, the sediment is eroded.
                p.resus(i) = resus_rate * ((tau - tau_crit)/tau_crit)^alpha; % resus flux (same unit as the resus rate, difference not entirely clear)
            end
            b_idx = b_idx + p.Zn-1 -i;
        end
    end
    
    
    
    p.resus_benth = p.resus;
end


% calculating light intensity in the pelagic
for j = 1:p.Xn-1 % pooling over each column (j = column nr)
    int = 0;
    for i = 1:p.Zn - temp % traversing down one column (i = row nr)
        
        dz_middle = p.Z_vol(i,j) - p.Z(i,j); % distance from bottom of grid element (i-1,j) to the center of grid element (i,j)
        
        if(i==1) % width of grid element (i,j) in the z-direction
            dz = p.Z(i+1,j);
        else
            dz = p.Z(i+1,j) - p.Z(i,j);
        end
        % note that detritus is measured in [mgP/m^3] while the light
        % atteniation has the units [m^2 mg C^-1], this a division by the
        % light attenuation coefficient p.q is neccesary
        int_middle = int + (1/p.q.*p.kD.*D(index) +p.kA*A(index)+ p.kbg)*dz_middle; % this is the integral down to the center of grid element (i,j).
        int = int + (1/p.q.*p.kD.*D(index) +p.kA*A(index)+ p.kbg)*dz; % here we integrate to the bottom of grid element (i,j)
        I(index) =  p.I0 .* exp( -int_middle);
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

% Algae consume and respire nutrients, nutrients from detritus remineralize
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
        dRsdt(i) =  dRsdt(i) + (1-p.benth_recycling).*p.q.*p.lbg_benth*B(i);
        
        if(p.model_version == 6)
            % Algae sinks into the sediment and becomes benthic algae
            dAdt(b_idx) = dAdt(b_idx) - p.vA.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.volumes_cyl(b_idx);
            dBdt(i) = dBdt(i)     + p.vA.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.Area_bottom_cyl(i);
        end
        
        if(p.model_version == 5)
            % Algae sinks into the sediment, dies and is turned into
            % sedimented nutrients.
            dAdt(b_idx) = dAdt(b_idx) - p.vA.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.volumes_cyl(b_idx);
            dRsdt(i) = dRsdt(i)  + p.q*p.vA.*A(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.Area_bottom_cyl(i);
        end
        
        
        
        % Detritus sinks into the sediment,and becomes sedimented nutrients.
        dDdt(b_idx) = dDdt(b_idx) - p.vD.*D(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.volumes_cyl(b_idx);
        dRsdt(i) = dRsdt(i)    + p.vD.*D(b_idx)*2*pi*0.5*(p.X(end,i+1)^2 - p.X(end,i)^2)./p.Area_bottom_cyl(i);
        
        if(p.model_version == 6)
            %benthic algae is resuspended into the water as pelagic algae
            dBdt(i) =  dBdt(i) - p.resus_benth(i).*B(i);
            dAdt(b_idx) = dAdt(b_idx) +  p.resus_benth(i).*B(i).*p.Area_bottom_cyl(i)/p.volumes_cyl(b_idx);
        end
        
        
        % Detritus is resuspended into the water
        dRsdt(i) =  dRsdt(i) - p.resus(i).*Rs(i);
        dDdt(b_idx) = dDdt(b_idx) +  p.resus(i).*Rs(i).*p.Area_bottom_cyl(i)/p.volumes_cyl(b_idx);
        
        
        %%%%% net flux of dissolved nutrients into the lake from the sediment and consumption by the benthic algae   %%%%%%%%%%%%%%%%%%%%%%%
        
        nutrient_limited_flux =  p.q.*B(i)'.*p.Gmax_benth.*(Rd(b_idx)./(Rd(b_idx)+p.M_benth)); % [mgP/ (m^2 day)]
        light_limited_flux = p.q.*p.Gmax./p.kB*log((p.H_benth + I(b_idx))./(p.H_benth +I(b_idx).*exp(-p.kB.*B(i)'))); % [mgP/ (m^2 day)]
        
        benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,
        net_flux = p.r.*Rs(i) +  p.benth_recycling.*p.q.*p.lbg_benth*B(i)' - benth_P_consumption;
        
        dRdt(b_idx) =   dRdt(b_idx) + net_flux.*p.Area_bottom_cyl(i)./p.volumes_cyl(b_idx); % the net flux here has the units [mg P/ (m^2 day)], and is multiplied by the area of the bttom segment to yield
        % a change [mg P/ day]. Division by the element volume yields the sought change in concentration.
        
        b_idx = b_idx + p.Zn-1 -i;
    end
end

%% Output

output = [dAdt(:); dRdt(:); dDdt(:); dRsdt(:); dBdt(:)]; % The output is a column vector with the derivatives of the state variables at each grid point.


end
