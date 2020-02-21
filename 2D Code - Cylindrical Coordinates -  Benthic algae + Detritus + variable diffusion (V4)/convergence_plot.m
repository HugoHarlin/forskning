
% Loads data and plots the convergence rate of the numerical method.

%% loading data
res_4 = load("2D_convergence_res_V4_dx_Xn_5_Zn_5");
res_16 = load("2D_convergence_res_V4_dx_Xn_17_Zn_17");
res_32 = load("2D_convergence_res_V4_dx_Xn_33_Zn_33");

%% Calculating the phytoplankton concentration at the point (res_4.p.X_vol(1,1), res_4.p.Z_vol(1,1))
% for three different resolutions. Because the points don't overlap whe
% nthe grid is refined, linear interpolation is used to calculated the
% concentration at the sought point.
val_res_4 = res_4.Y_t(end,4^2+1);
val_res_16 = 0.25*(res_16.Y_t(end,16^2+1) + res_16.Y_t(end,16^2 + 2) + res_16.Y_t(end,16^2 + 16) + res_16.Y_t(end,16^2 + 17));
val_res_32 = 0.25*(res_32.Y_t(end,32^2 + 34) + res_32.Y_t(end,32^2 + 35) + res_32.Y_t(end,32^2 + 32+24) + res_32.Y_t(end,32^2 + 32+35));

%% Extracting detritus
D_res_4   = res_4.Y_t(end,2*(res_4.p.Xn-1)*(res_4.p.Zn-1)+1 : 3*(res_4.p.Xn-1)*(res_4.p.Zn-1));
D_res_16  = res_16.Y_t(end,2*(res_16.p.Xn-1)*(res_16.p.Zn-1)+1 : 3*(res_16.p.Xn-1)*(res_16.p.Zn-1));
D_res_32  = res_32.Y_t(end,2*(res_32.p.Xn-1)*(res_32.p.Zn-1)+1 : 3*(res_32.p.Xn-1)*(res_32.p.Zn-1));

D_res_4  = reshape(D_res_4, [res_4.p.Xn-1, res_4.p.Zn-1]);
D_res_16  = reshape(D_res_16, [res_16.p.Xn-1, res_16.p.Zn-1]);
D_res_32  = reshape(D_res_32, [res_32.p.Xn-1, res_32.p.Zn-1]);

D_res_4 = D_res_4';
D_res_16 = D_res_16';
D_res_32 = D_res_32';

%% calculation of relative error, based on the norm of the state variabless.

p_top_left = log(abs((val_res_32-val_res_16)/((val_res_16- val_res_4))))/log(2);



%% plotting results
plot(rel_vec);

