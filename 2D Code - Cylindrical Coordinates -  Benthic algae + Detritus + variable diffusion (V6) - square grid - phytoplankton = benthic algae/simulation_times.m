figure(1)
% efficient code
resolution = [2,10,15,20,25] ;
simtime_eff = [0.186181686369310,1.546511233551906 ,9.071942577220298, 38.323501643031640, 1.135326657153579e+02];
simtime_slow = [0.254907852746731, 7.653903617805788, 32.079728333625680, 1.001394912163175e+02, 2.385287524154990e+02];
% the resolution of the slow sime time is 1 lower than it should be...

plot(resolution, simtime_eff,resolution, simtime_slow);
xlabel("resolution [nr gridpoints in x and z dimension]");
ylabel("simulation time [seconds]");


x = 40;
T = 0.4*x^2 -6.2*x +14;

T = T/(60*60)



%% simulation time for tolerances set to 10e-6

resolution = [ 5, 10, 15, 20, 25, 30];
simtime = [0.281986968795345, 1.139946692184012, 6.403818193288787, 16.999090506956673, 1.509482438652069e+02, 939];
ordo_n3 = resolution.^3;
ordo_n3 = ordo_n3*max(900/ordo_n3(end));
plot(resolution, simtime, resolution, ordo_n3)