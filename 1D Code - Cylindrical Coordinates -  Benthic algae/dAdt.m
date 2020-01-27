function [output] = dAdt(t,Y,p,I)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface laZer and so forth.


%% reformatting input and calculation of light intensity etc.
% separating the state variables from Y
A = Y(1:(p.Zn-1));
Rd = Y(p.Zn:2*(p.Zn-1));
Rs = Y(2*(p.Zn-1)+1);
B = Y(2*(p.Zn-1)+2);

% Calculation of the light intensity at each grid element center
I = zeros(p.Zn-1,1);
integ = 0;
for i=1:p.Zn-1
    integ = integ + p.Z(i) *A(i);
    I(i) = p.I0 * exp(-p.k*integ -p.kbg*p.Z(i));
end

G = p.Gmax .* min(Rd./(Rd+p.M), I./(I+p.H)); % New source term. The production is not colimited by light
%  and nutrients as in the original case, and is instead only limited by
%  the most limiting resource. This design choice was made because the
%  growth in the benthic layer is designed this way, and keeping the colimited phytoplankton
% would imply a discrepancy in the underlying assumptions of growth which is
%  inconsistent.

%% Diffusion and Convection

% diffusion of dissolved Nutrients
dRdt = p.diffusion_matrix*Rd;

% diffusion and sinking of algae
dAdt = p.dynam_matrix*A;

% the sinking algae that reaches the bottom is turned into sedimented
% nutrients.
dRsdt = A(end)*p.v*p.q; %

%% Source terms
% Algae grow with net rate G - p.lbg, [mgC/(m^3 day)]
dAdt = dAdt + (G - p.lbg).*A;

% correspoding decrease in nutrients [mgP/(m^3 day)]
dRdt = dRdt - p.q.*(G-p.lbg).*A;

% Benthic Algae respire/die and a portion p.benth_recycling is recycled for
% the immediate use of the benthic algae. The rest is introduced into the
% lake as dissolved nutrients.
dRsdt = dRsdt + (1-p.benth_recycling).*p.lbg_benth.*p.q_benth.*B'; 

% The nutrients in the sediment remineralize back into the lake. The
% benthic algae gobble up as much as they can, the rest is introduced into
% the water as dissolved nutrients.
dRsdt = dRsdt - p.r*Rs;

% Growth of bentic algae, which is limitied by either light or
% nutrients via a min-function. The Benthic algae respire nitrogen at
% a constant rate lbg_benth. A part of those nutrients, benth_recycling, can be
% immediately reused by the benthic algae for growth. The benthic algae can
% also use the dissolved nutrients in the water closest to the sediment for
% growth.

nutrient_limited_flux = Rs*p.r + p.q_benth*B*p.Gmax_benth*(Rd(end)/(Rd(end)+p.M_benth)) ...
     + p.q_benth*p.benth_recycling*p.lbg_benth*B; % [mgP/(m^2 day)]

light_limited_flux = p.q_benth.*p.Gmax_benth./p.kB.*log((p.H + I(end))/(p.H + I(end).*exp(-p.k*B))); % [mgP/ (m^2 day)]
benth_P_consumption = min(nutrient_limited_flux, light_limited_flux); %This is the total consumption of dissolved nutrients,

% Net flux pertains to the total amount of nutrients released from the
% sediemt and the benthic algae minus the ammount that the benthic
% algae consumes. The net flux is thus the total flux of nutrients back
% into the lake from the benthic algae and sediment layer.
net_flux = p.benth_recycling*B*p.lbg_benth*p.q_benth + p.r*Rs(end) - benth_P_consumption;

% net flux of mineralized nutrients in the water at the bottom.
dRdt(end) = dRdt(end) + net_flux/p.deltaZ; 

% Growth of benthic algae [mgC/(m^2 day)]
dBdt = benth_P_consumption*(1/p.q_benth) - p.lbg_benth*B; % The 1/p.q_benth factor is neccesary since benth_P_consumption has the units [mgP/(m^2 day)] 
                                          % and dBdt is measured in
                                          % [mgC/(m^2 day)].
%% reshaping matrices for output
% transpose of matrices in order to use colon notation to reshape to vector form.
dAdt = dAdt';
dRdt = dRdt';
output = [dAdt(:); dRdt(:); dRsdt; dBdt]; % The output is a column vector with the derivatives of the state variables at each grid point.
end







