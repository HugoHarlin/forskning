% simtime plot for v4 for increasing resolution.

sim_t = [7.452164380844091, 16.104100258442788, 29.376547164988022, 51.617634134431825, 81.803509453929760, 1.585813861387583e+02, 2.353180965015956e+02, 3.631274589352391e+02, 5.327770678431460e+02 ];
res   = [10, 12, 14, 16, 18, 20, 22, 24, 26 ];

figure(11)
plot(res,sim_t);

x = linspace(0,100,100);
y = 0.14.*x.^3 -4.5.*x.^2 + 54.*x -2.2e+02;

figure(12)
plot(x,y);

