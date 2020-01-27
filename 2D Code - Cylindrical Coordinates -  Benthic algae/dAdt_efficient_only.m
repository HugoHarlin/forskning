function [output] = dAdt_efficient_only(t,Y,p)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface layer and so forth.

% function is used as input to ode15s

% ONLY CYLINDRICAL COORDINATES - no cartesian implementation in this script


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

% Bottom sediment interaction and benthic algal growth
i = p.Zn-1;
for j = 1:p.Xn-1 % xi
    % Bethic algae net growth
    nutrient_limited_growth = 1/p.q_benth*Rs(j)*p.r + p.Gmax_benth*(Rd(i,j)/(Rd(i,j)+p.M_benth))*B(j) + ...
        p.benth_recycling*p.lbg_benth*B(j);
    
    light_limited_growth = p.Gmax./p.kB.*log((p.H_benth + I(i,j))/(p.H_benth + I(i,j).*exp(-p.k.*B(j))));
    benthic_growth = min(nutrient_limited_growth, light_limited_growth);
    dBdt(j) = dBdt(j) + benthic_growth - p.lbg_benth*B(j);
    
    %     if(light_limited_growth < nutrient_limited_growth)
    %         disp("light limited growth");
    %     else
    %         disp("nutrient limited growth");
    %     end
    
    % Nutrients in the sediment remineralize. The benthic algae gobble up
    % as much as they can, the rest is introduced into the lake as
    % dissolved nutrients.
    dRsdt(j) =  dRsdt(j) - 2*pi*p.W/(p.Xn-1).*(j-0.5).*Rs(j)*p.r*sqrt(p.dX_dXi_preCalc(i,j,4)^2 + p.dZ_dXi_preCalc(i,j,4)^2);
    
    % Benthic Algae respire/die and a portion p.benth_recycling is bound in
    % particulate form and is introduced into the sediment. This is handled
    % after the division by the element area below.
    
    % Algae sink into the sediment, dies and instantly becomes sedimented
    % nutrients.
    dRsdt(j) =  dRsdt(j) + p.death_rate*p.q*2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,4)*A(i,j);
    
    
        %     %%%%% net flux of dissolved nutrients into the lake from the sediment and consumption by the benthic algae   %%%%%%%%%%%%%%%%%%%%%%%
        %net_flux =   p.benth_recycling*B(j)*p.lbg_benth*p.q_benth + p.r*Rs(j) - p.q*benthic_growth;
        
        nutrient_limited_flux = Rs(j)*p.r + p.q_benth*p.Gmax_benth*(Rd(i,j)/(Rd(i,j)+p.M_benth))* ...
            B(j) + p.q_benth*p.benth_recycling*p.lbg_benth*B(j); % [mgP/ (m^2 day)]
        
        light_limited_flux = p.q_benth.*p.Gmax./p.kB*log((p.H_benth + p.I(i,j))/(p.H_benth + p.I(i,j).*exp(-p.k.*B(j)))); % [mgP/ (m^2 day)]
        benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,
        
        net_flux = p.benth_recycling*B(j)*p.lbg_benth*p.q_benth + p.r*Rs(j) - benth_P_consumption;
        
        dRdt(i,j) =   dRdt(i,j) + 2*pi*p.W/(p.Xn-1).*(j-0.5).*net_flux*sqrt(p.dX_dXi_preCalc(i,j,4)^2 + p.dZ_dXi_preCalc(i,j,4)^2);
    
    
end

%% Division by element area
% This is only done for the diffusion and advection contributions.
% the volume is a factor in the source and sink terms and thus
% cancels when the l.h.s. and r.h.s are divided by the volume V.
% This cancellation does not occur for the convection and diffusion terms
% since we utiliyed the divergence theorem to rewrite the integrals as
% line integrals, and the volume of the correspoding element does not
% appear as a factor in that case.

dAdt = dAdt./p.volumes_cyl;
dRdt = dRdt./p.volumes_cyl;
dRsdt = dRsdt./p.L_bottom_cyl;

%% Source terms
% Algae grow with net rate g - p.lbg, note that algae is measured in
% [mgC/m^3]
dAdt = dAdt + (G - p.lbg).*A;

% correspoding decrease in nutrients [mgP/m^3]
dRdt = dRdt - p.q.*(G-p.lbg).*A;

% Benthic Algae respire/die and a portion p.benth_recycling is bound in
% particulate form and is introduced into the sediment.
dRsdt = dRsdt + (1-p.benth_recycling).*p.lbg_benth.*p.q_benth.*B';

%% reshaping matrices for output
% transpose of matrices in order to use colon notation to reshape to vector form.
dAdt = dAdt';
dRdt = dRdt';
output = [dAdt(:); dRdt(:); dRsdt(:); dBdt(:)]; % The output is a column vector with the derivatives of the state variables at each grid point.

    % diffusion of algae and dissolved nutrients, sinking of algae
    %output(1:2*(p.Xn-1)*(p.Zn-1)) = output(1:2*(p.Xn-1)*(p.Zn-1)) + p.S(1:2*(p.Xn-1)*(p.Zn-1),1:2*(p.Xn-1)*(p.Zn-1))*Y(1:2*(p.Xn-1)*(p.Zn-1));
    output = output + p.S*Y;
end





