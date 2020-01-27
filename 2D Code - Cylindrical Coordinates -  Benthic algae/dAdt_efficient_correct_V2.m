function [output] = dAdt_efficient_correct_V2(t,Y,p)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface layer and so forth.

% function is used as input to ode15s

% V2 implies that the boundary conditions at the bottom are changed.
% The benthic algae can only access nutrients from the water column
% directly above the benthic layer, and cannot use the nutrients directly
% in the sediment or the remineralized sediment directly.

efficient = true;

%% reformatting input and calculation of light intensity etc.
% separating the state variables from y
A = Y(1:(p.Xn-1)*(p.Zn-1));
Rd = Y((p.Xn-1)*(p.Zn-1)+1 : 2*(p.Xn-1)*(p.Zn-1));
Rs = Y(2*(p.Xn-1)*(p.Zn-1) +1 : 2*(p.Xn-1)*(p.Zn-1) + (p.Xn-1));
B = Y(2*(p.Xn-1)*(p.Zn-1) + (p.Xn-1) +1 : end);

A = reshape(A,[(p.Xn-1) (p.Zn-1)]);
Rd = reshape(Rd,[(p.Xn-1) (p.Zn-1)]);
A = A';
Rd = Rd';

% Calculation of the light intensity at each grid element center
I = zeros(p.Zn-1,p.Xn-1);
for j = 1:p.Xn-1
    integral = 0;
    for i=1:p.Zn-1
        integral = integral + p.Z(i,j) *A(i,j);
        I(i,j) = p.I0 * exp(-p.k*integral -p.kbg*p.Z(i,j));
    end
end

p.I = I;
% Calculation of the algal production G(R,I)
%G =  p.Gmax .* Rd./(Rd+p.M) .* I./(I+p.H); This is the old source term.

G = p.Gmax .* min(Rd./(Rd+p.M), I./(I+p.H)); % New source term. The production is not colimited by light
%  and nutrients as in the original case, and is instead only limited by
%  the most limiting resource. This design choice was made because the
%  growth in the benthic layer is designed this way, and keeping the colimited phytoplankton
% would imply a discrepancy in the underlying assumptions of growth which is
%  inconsistent.

%% Diffusion and Convection, Bottom boundary condition
dAdt = zeros(p.Zn-1, p.Xn-1);
dRdt = zeros(p.Zn-1, p.Xn-1);
dRsdt = zeros(1,p.Xn-1);
dBdt = zeros(1,p.Xn-1);

if(~efficient)
    % Diffusion terms in cylindrical coordinates with grid transform. Boundary
    % conditions are incorporated in the integral functions defined at the end of the script.
    for j = 1:p.Xn-1 % xi
        for i = 1:p.Zn-1 % eta
            % Diffusion and sinking of algae with  BC:s
            
            dAdt(i,j) = integral_A_1(i,j,p,A) + integral_A_2(i,j,p,A) + ...
                integral_A_3(i,j,p,A) + integral_A_4(i,j,p,A);
            
            %diffusion of nutrients with BC:s
            dRdt(i,j) = integral_Rd_1(i,j,p,Rd) + integral_Rd_2(i,j,p,Rd) + ...
                integral_Rd_3(i,j,p,Rd) + integral_Rd_4(i,j,p,Rd,Rs,B);
            
        end
    end
end

% Bottom sediment interaction and benthic algal growth
i = p.Zn-1;
for j = 1:p.Xn-1 % xi
    % Bethic algae net growth
    
    nutrient_limited_growth =  p.Gmax_benth*B(j)*(Rd(i,j)/(Rd(i,j)+p.M_benth));
    light_limited_growth = p.Gmax_benth./p.kB.*log((p.H_benth + I(i,j))/(p.H_benth + I(i,j).*exp(-p.k.*B(j))));
    
    dBdt(j) = dBdt(j) + min(nutrient_limited_growth, light_limited_growth) - p.lbg_benth*B(j);
    
    
    % Nutrients in the sediment remineralize into the lake.
     dRsdt(j) =  dRsdt(j) -Rs(j)*p.r;
    
    % Benthic Algae respire/die and a portion p.benth_recycling is bound in
    % particulate form and is introduced into the sediment. This is handled
    % after the division by the element area below.
    
    % Algae sink into the sediment, dies and instantly becomes sedimented
    % nutrients.
    dRsdt(j) =  dRsdt(j) + p.v*p.death_rate*A(i,j)*p.q/p.L_bottom_cyl(j)*2*pi*p.W/(p.Xn-1).*(j-0.5).*p.dX_dXi_preCalc(i,j,4);
    
    
    if(efficient)
        %     %%%%% net flux of dissolved nutrients into the lake from the sediment and consumption by the benthic algae   %%%%%%%%%%%%%%%%%%%%%%%
        %net_flux =   p.benth_recycling*B(j)*p.lbg_benth*p.q_benth + p.r*Rs(j) - p.q*benthic_growth;
        
        nutrient_limited_flux =  p.q_benth*B(j)*p.Gmax_benth*(Rd(i,j)/(Rd(i,j)+p.M_benth)); % [mgP/ (m^2 day)]
        
        light_limited_flux = p.q_benth.*p.Gmax./p.kB*log((p.H_benth + p.I(i,j))/(p.H_benth + p.I(i,j).*exp(-p.k.*B(j)))); % [mgP/ (m^2 day)]
        benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,
        
        
        
        net_flux = p.r*Rs(j) + p.q_benth.*p.lbg_benth*B(j) - benth_P_consumption;
        
        dRdt(i,j) =   dRdt(i,j) + net_flux*p.L_bottom_cyl(j); % the net flux here has the units [mg P/ (m^2 day)], and is multiplied by the area of the bttom segment to yield
                                                              % a change [mg P/ day]. Division by the element volume yields the sought change in concentration.
                                    
                                                
                                                             
    end
    
    
end

%% Division by element volume
% This is only done for the diffusion and advection contributions.
% the volume is a factor in the source and sink terms and thus
% cancels when the l.h.s. and r.h.s are divided by the volume V.
% This cancellation does not occur for the convection and diffusion terms
% since we utiliyed the divergence theorem to rewrite the integrals as
% line integrals, and the volume of the correspoding element does not
% appear as a factor in that case.


dAdt = dAdt./p.volumes_cyl;
dRdt = dRdt./p.volumes_cyl;
% dRsdt = dRsdt./p.L_bottom_cyl;

%% Source terms
% Algae grow with net rate g - p.lbg, note that algae is measured in
% [mgC/m^3]
dAdt = dAdt + (G - p.lbg).*A;

% correspoding decrease in nutrients [mgP/m^3]
dRdt = dRdt - p.q.*(G-p.lbg).*A;

% Benthic Algae respire/die and a portion p.benth_recycling is bound in
% particulate form and is introduced into the sediment.
%dRsdt = dRsdt + (1-p.benth_recycling).*p.lbg_benth.*p.q_benth.*B';

%% reshaping matrices for output
% transpose of matrices in order to use colon notation to reshape to vector form.
dAdt = dAdt';
dRdt = dRdt';
output = [dAdt(:); dRdt(:); dRsdt(:); dBdt(:)]; % The output is a column vector with the derivatives of the state variables at each grid point.

if(efficient)
    % diffusion of algae and dissolved nutrients, sinking of algae
    output = output + p.S*Y;
end

end



%% Support functions
%%% Line integrals Algae %%%
function [output] = integral_A_1(i,j,p,A) % Integral of right side.  j=xi, i=eta

if(j ~= p.Xn-1) % we are not on the right border.
    dA_dXi = A(i,j+1) - A(i,j);
    output = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi - ...
        p.dZ_dXi_preCalc(i,j,1) * dA_dEta_1(i,j,p,A));
else
    output = 0; % BC:s imply that there is no flux on the right border.
end
end
function [output] = integral_A_2(i,j,p,A) % Integral of top side.    j=xi, i=eta
if(i ~=1)
    dA_dEta = A(i,j) - A(i-1,j);
    
    % Diffusion
    output = 2*pi*p.W/(p.Xn-1).*(j-0.5).*(p.dx * p.dZ_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)* ...
        dA_dXi_2(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dZ_dEta_preCalc(i,j,2))* ...
        dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,2)^2 + p.dz*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Sinking
    output = output + 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,2).*0.5.*(A(i,j)+A(i-1,j));
else
    output = 0; % BC:s imply that there is no flux on the top border.
end
end
function [output] = integral_A_3(i,j,p,A) % Integral of left side.   j=xi, i=eta
if(j ~= 1) % we are not on the left border.
    dA_dXi = A(i,j)-A(i,j-1);
    output = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/p.dX_dXi_preCalc(i,j,3) * ...
        (p.dZ_dEta_preCalc(i,j,3) .* dA_dXi - ...
        p.dZ_dXi_preCalc(i,j,3) * dA_dEta_3(i,j,p,A));
else
    output = 0; % BC:s imply that there is no flux on the left border.
end
end
function [output] = integral_A_4(i,j,p,A) % Integral of bottom side. j=xi, i=eta
if(i ~=p.Zn-1) % No diffusive flux through the bottom border.
    dA_dEta = A(i+1,j) - A(i,j);
    
    %Diffusion
    output = -2*pi*p.W/(p.Xn-1).*(j-0.5) .* (p.dx * p.dZ_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        dA_dXi_4(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,4)*p.dZ_dEta_preCalc(i,j,4))* ...
        dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,4)^2 + p.dz*p.dX_dXi_preCalc(i,j,4)^2));
    
    % Sinking
    output = output - 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,4).*0.5.*(A(i+1,j)+A(i,j));
else
    % algae that sink into the bottom become sedimented algae
    output = - p.death_rate*2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,4)*A(i,j);
end
end

%%% Differential Quadratures Algae %%%
function [dA_dEta] = dA_dEta_1(i,j,p,A) % j=xi, i=eta. dA_dEta for integral 1
if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
    dA_dEta = 0.25*(A(i,j) + A(i,j+1) - A(i-1,j) - A(i-1,j+1));
elseif(i ==  1) %  top border
    dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i,j) - A(i,j+1));
else % in between top and bottom border. we can use a standard quadrature
    dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i-1,j) - A(i-1,j+1));
end
end
function [dA_dEta] = dA_dEta_3(i,j,p,A) % j=xi, i=eta. dA_dEta for integral 3
if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
    dA_dEta = 0.25*(A(i,j) + A(i,j-1) - A(i-1,j) - A(i-1,j-1));
elseif(i ==  1) %  top border
    dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i,j) - A(i,j-1));
else % in between top and bottom border. we can use a standard quadrature
    dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i-1,j) - A(i-1,j-1));
end
end
function [dA_dXi] = dA_dXi_2(i,j,p,A)   % j=xi, i=eta. dA_dXi  for integral 2
if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
    dA_dXi = 0.25*(A(i,j) + A(i-1,j) - A(i,j-1) - A(i-1,j-1));
elseif(j ==  1) %  left
    dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j) - A(i-1,j));
else % in between left and right border. we can use a standard quadrature
    dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j-1) - A(i-1,j-1));
end
end
function [dA_dXi] = dA_dXi_4(i,j,p,A)   % j=xi, i=eta. dA_dXi  for integral 4
if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
    dA_dXi = 0.25*(A(i,j) + A(i+1,j) - A(i,j-1) - A(i+1,j-1));
elseif(j ==  1) %  left
    dA_dXi = 0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j) - A(i+1,j));
else % in between left and right border. we can use a standard quadrature
    dA_dXi = 0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j-1) - A(i+1,j-1));
end
end

%%% Line integrals Dissolved Phosphorus %%%
function [output] = integral_Rd_1(i,j,p,Rd)    % Integral of right side.  j=xi, i=eta
if(j ~= p.Xn-1) % we are not on the right border.
    dRd_dXi = Rd(i,j+1) - Rd(i,j);
    output = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        (p.dZ_dEta_preCalc(i,j,1) .* dRd_dXi - ...
        p.dZ_dXi_preCalc(i,j,1) * dRd_dEta_1(i,j,p,Rd));
else
    output = 0; % BC:s imply that there is no flux on the right border.
end
end
function [output] = integral_Rd_2(i,j,p,Rd)    % Integral of top side.    j=xi, i=eta
if(i ~=1)
    dRd_dEta = Rd(i,j) - Rd(i-1,j);
    
    % Diffusion
    output = 2*pi*p.W/(p.Xn-1).*(j-0.5).*(p.dx * p.dZ_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)* ...
        dRd_dXi_2(i,j,p,Rd) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dZ_dEta_preCalc(i,j,2))* ...
        dRd_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,2)^2 + p.dz*p.dX_dXi_preCalc(i,j,2)^2));
else
    output = 0; % BC:s imply that there is no flux on the top border.
end
end
function [output] = integral_Rd_3(i,j,p,Rd)    % Integral of left side.   j=xi, i=eta
if(j ~= 1) % we are not on the left border.
    dRd_dXi = Rd(i,j)-Rd(i,j-1);
    output = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/p.dX_dXi_preCalc(i,j,3) * ...
        (p.dZ_dEta_preCalc(i,j,3) .* dRd_dXi - ...
        p.dZ_dXi_preCalc(i,j,3) * dRd_dEta_3(i,j,p,Rd));
else
    output = 0; % BC:s imply that there is no flux on the left border.
end
end
function [output] = integral_Rd_4(i,j,p,Rd,Rs,B) % Integral of bottom side. j=xi, i=eta
if(i ~=p.Zn-1)
    dRd_dEta = Rd(i+1,j) - Rd(i,j);
    
    %Diffusion
    output = -2*pi*p.W/(p.Xn-1).*(j-0.5) .* (p.dx * p.dZ_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        dRd_dXi_4(i,j,p,Rd) - 1/( p.dX_dXi_preCalc(i,j,4)*p.dZ_dEta_preCalc(i,j,4))* ...
        dRd_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,4)^2 + p.dz*p.dX_dXi_preCalc(i,j,4)^2));
else
    % Here the boundary condition applies, which is handled separately
    nutrient_limited_flux = Rs(j)*p.r + p.q_benth*p.Gmax_benth*(Rd(i,j)/(Rd(i,j)+p.M_benth))* ...
        B(j) + p.q_benth*p.benth_recycling*p.lbg_benth*B(j); % [mgP/ (m^2 day)]
    
    light_limited_flux = p.q_benth.*p.Gmax./p.kB*log((p.H_benth + p.I(i,j))/(p.H_benth + p.I(i,j).*exp(-p.k.*B(j)))); % [mgP/ (m^2 day)]
    benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,
    
    net_flux = p.benth_recycling*B(j)*p.lbg_benth*p.q_benth + p.r*Rs(j) - benth_P_consumption;
    output = 2*pi*p.W/(p.Xn-1).*(j-0.5).*net_flux*sqrt(p.dX_dXi_preCalc(i,j,4)^2 + p.dZ_dXi_preCalc(i,j,4)^2);
end
end

%%% Differential Quadratures Phosphorus %%%
function [dRd_dEta] = dRd_dEta_1(i,j,p,Rd) % j=xi, i=eta. dA_dEta for integral 1
if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
    dRd_dEta = 0.25*(Rd(i,j) + Rd(i,j+1) - Rd(i-1,j) - Rd(i-1,j+1));
elseif(i ==  1) %  top border
    dRd_dEta = 0.25*(Rd(i+1,j) + Rd(i+1,j+1) - Rd(i,j) - Rd(i,j+1));
else % in between top and bottom border. we can use a standard quadrature
    dRd_dEta = 0.25*(Rd(i+1,j) + Rd(i+1,j+1) - Rd(i-1,j) - Rd(i-1,j+1));
end
end
function [dRd_dEta] = dRd_dEta_3(i,j,p,Rd) % j=xi, i=eta. dA_dEta for integral 3
if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
    dRd_dEta = 0.25*(Rd(i,j) + Rd(i,j-1) - Rd(i-1,j) - Rd(i-1,j-1));
elseif(i ==  1) %  top border
    dRd_dEta = 0.25*(Rd(i+1,j) + Rd(i+1,j-1) - Rd(i,j) - Rd(i,j-1));
else % in between top and bottom border. we can use a standard quadrature
    dRd_dEta = 0.25*(Rd(i+1,j) + Rd(i+1,j-1) - Rd(i-1,j) - Rd(i-1,j-1));
end
end
function [dRd_dXi] = dRd_dXi_2(i,j,p,Rd)   % j=xi, i=eta. dA_dXi  for integral 2
if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
    dRd_dXi = 0.25*(Rd(i,j) + Rd(i-1,j) - Rd(i,j-1) - Rd(i-1,j-1));
elseif(j ==  1) %  left
    dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i-1,j+1) - Rd(i,j) - Rd(i-1,j));
else % in between left and right border. we can use a standard quadrature
    dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i-1,j+1) - Rd(i,j-1) - Rd(i-1,j-1));
end
end
function [dRd_dXi] = dRd_dXi_4(i,j,p,Rd)   % j=xi, i=eta. dA_dXi  for integral 4
if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
    dRd_dXi = 0.25*(Rd(i,j) + Rd(i+1,j) - Rd(i,j-1) - Rd(i+1,j-1));
elseif(j ==  1) %  left
    dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i+1,j+1) - Rd(i,j) - Rd(i+1,j));
else % in between left and right border. we can use a standard quadrature
    dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i+1,j+1) - Rd(i,j-1) - Rd(i+1,j-1));
end
end





