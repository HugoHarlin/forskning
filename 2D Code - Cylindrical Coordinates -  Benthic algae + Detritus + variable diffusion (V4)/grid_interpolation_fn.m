function [X_vol_new,Y_vol_new] = grid_interpolation_fn(p,res_x,res_z)
%% Mesh generator
%Generates finer mesh for interpolation of results, makes for nicer plots.
%
% Hugo Harlin 2019

% Quantities relating to system size
p.Xn_new = res_x+1; % Number of grid-points (width)
p.Yn_new = res_z+1; % Number of grid-points (depth)

% Lake Mesh, with an increasing depth from Lmin at the shore to Lmax
% at the center of the lake ( slope = alpha* (Lmin - Lmax)/W  ).
% (0,0) is placed at the center of the lake at the surface, y-dim is facing
% downward and x-dim is facing towards the lake edge

p.X_new = zeros(p.Yn_new, p.Xn_new); % Mesh-spacing in x-dimension
p.Y_new = zeros(p.Yn_new, p.Xn_new); % Mesh-spacing in y-dimension

for i=1:1:p.Yn_new
    p.X_new(i,:) = [0:p.W/(p.Xn_new-1):p.W]; % even spacing of the grid in x-dimension
end

% the grid is compressed in y-dimension, with depth Lmin at the
% shore and Lmax at the center of the lake.
for i=1:p.Yn_new
    for j = 1:p.Xn_new
        p.Y_new(i,j) = (p.Lmax/(p.Yn_new-1))*(i-1)*(1 + (p.Lmin/p.Lmax -1)* (p.X_new(i,j)/p.W).^(p.alpha));
    end
end

% coordinates of the center of each mesh quadrilateral
% returns the coordinates of the center of each volume element in the mesh
X_vol_new = zeros(p.Yn_new-1, p.Xn_new-1);
Y_vol_new = zeros(p.Yn_new-1, p.Xn_new-1);


for j=1:p.Xn_new-1
    for i =1:p.Yn_new-1
        X_vol_new(i,j) = (p.X_new(i,j) + p.X_new(i+1,j) + p.X_new(i,j+1) + p.X_new(i+1,j+1) ) /4;
        Y_vol_new(i,j) = ( p.Y_new(i,j) +  p.Y_new(i+1,j) +  p.Y_new(i,j+1) +  p.Y_new(i+1,j+1) ) /4;
    end
end

end

