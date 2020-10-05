% håkansson resuspension 2008 & 2003

% This script calculates the resuspension from erosion-transport (ET) areas
% according to the formulations in articles
% " A Dynamic Mass-balance Model for Phosphorus in Lakes with a Focus on Criteria for Applicability and Boundary Conditions" (2008)
% " A Dynamic Model to Predict Suspended Particulate Matter in Lakes" (2003)
% by Håkansson.

% The formula below is a combination of the very similar approaches in the
% articles above.

% in the 2008 model the average wind speed is assumed to be the ref value,
% and the sediment age Tet is calculated.

% In the 2003 paper this assumption is not used and the surface wind speed
% is used in the calculation, but the sediment age Tet is instead assumed to be 12 months.

% The smallest value the resuspension can obtain is 0.27%/day. (Tet = 12,
% wind speed = ref value = 3.27m/s


%%%%%%%%%%%%%%%%%%%%%%% The formula  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total resuspension from the ET sediment, given by the equations in table
% 4 under the ET-areas section in the 2008 paper (FETSW and FETDW). Res_tot = (1-Vd)Rres + Vd*Rres = Rres [1/month]
% Rres is found further down in table 4 and is given by

% Rres = 1/Tet,  [1/month]  ( using the wind speed modelling from 2003 paper Rres = (wind_speed/ref_wind_speed)^2/Tet )

% where Tet is the sediment age in months and is given by equation 9 (in
% 2003 paper it is assumed to be 12 months, i.e. dr is assumed to be 0.26)

% Tet = 12(dr/0.26) if dr < 0.26
% Tet = 12(0.26/dr) if de > 0.26

% where dr is the dynamic ratio and is defined in table 4

% dr = sqrt((lake area)*10^-6)/dm, where dm is the mean lake depth
% dr is dimensionless
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dm = 15; % mean lake depth
radius = 2000; % lake radius in m
area = pi*radius^2; % lake surface are in m^2

%%%% data från Vänern (wikipedia) %%%%
%area = 5650*10^6; % vänerns yta [m^2]
%dm = 27; % vänerns medeldjup [m]
%%%%%%%%%%%%%%%%

wind = 3.27; % average surface wind speed in m/s
ref_wind_speed = 3.27; % reference wind speed (where does Håkansson get this? /H)


dr = sqrt(area*10^-6)/dm; % sediment age in months
disp("dynamic ratio: " + num2str(dr) + " [1]");


if(dr < 0.26)
    tet = 12*dr/0.26;
else
    tet = 12*0.26/dr;
end

resus = (wind/ref_wind_speed)^2/tet; % resuspension. unit [1/month]
resus = resus/31; % change of units to [1/day]

disp("resuspension: " + num2str(resus*100) + " [%/day]");


