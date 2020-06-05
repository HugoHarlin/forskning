% calculates and plots maximum benthic biomass as a function of available
% light



Gmax = 1.04; % max specific growth rate of benthic algae [day^-1]
H = 1.5;  % Half-saturation constant of light-dependent benthic algal production [micro-mol photons m^-2 s^-1]
kbgB = 1e-3; % light attenuation coefficient of benthic algae [m^2 mg C^-1]
lB = 0.1; % total specific loss of benthic alge

I = 1:10:1000; % light intensity at the surface of the benthic layer [

syms B
eqn = Gmax./kbgB.*log((H +I)./(H + I.*exp(-kbgB.*B))) - lB.*B;

max_B = zeros(length(I),1);

for i=1:length(I)
    i
    eqn = Gmax./kbgB.*log((H +I(i))./(H + I(i).*exp(-kbgB.*B))) - lB.*B;
    max_B(i) = vpasolve(eqn == 0,B,1e12);
end


loglog(I,max_B);
grid on
title("light limited max benthic biomass (10% loss rate)");
xlabel("incoming light intensity [micro-mol photons m^-2 s^-1]");
ylabel("Benthic biomass [mgC/m^2]");


