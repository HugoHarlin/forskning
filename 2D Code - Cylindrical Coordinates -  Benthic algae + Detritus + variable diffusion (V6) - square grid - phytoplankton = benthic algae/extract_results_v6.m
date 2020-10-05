function [r,A,Rd,D,Rs,B,ntot_0,n_algae_end,n_dissolved_end,n_detritus_end,n_sediment_end,n_benthic_end,ntot_end] = extract_results_v6(t,Y_t,p)
%% Extracts results from simulation

% unpackages and reformats the results from the input simulation,
% and returns the extracted and calculated values. 

%% Extracting results
ntot_0 = p.ntot_0;

% Time of plotted results
%time = length(t); % winter
%time = length(t)-183; % summer
%time = length(t)-274; % spring
%time = length(t)-91; % fall
time = length(t);

A  = Y_t(time,1:sum(1:p.Xn-1)); % phytoplankton
Rd = Y_t(time,sum(1:p.Xn-1)+1:2*sum(1:p.Xn-1)); % dissolved nutrients
D  = Y_t(time,2*sum(1:p.Xn-1) +1 : 3*sum(1:p.Xn-1)); % ditritus
Rs = Y_t(time,3*sum(1:p.Xn-1)+1 : 3*sum(1:p.Xn-1)+ p.Xn-1); % sedimented nutrients
B  = Y_t(time, 3*sum(1:p.Xn-1)+ p.Xn :  3*sum(1:p.Xn-1)+ 2*(p.Xn-1)); % benthic algae

r = p.resus(1);
%Calculating nutrient content at t = Tend
n_algae_end     = p.volumes_cyl.*A'*p.q;
n_dissolved_end = p.volumes_cyl.*Rd';
n_detritus_end  = p.volumes_cyl.*D';
n_sediment_end  = Rs.*p.Area_bottom_cyl;
n_benthic_end   = B.*p.q.*p.Area_bottom_cyl;
ntot_end        = sum(sum(n_algae_end)) + sum(sum(n_dissolved_end)) + sum(sum(n_detritus_end)) + sum(n_sediment_end) + sum(n_benthic_end);

% creation of refined grid for interpolation.
%[X_vol_new,Y_vol_new] = grid_interpolation_fn(p,100,100);
end

