function [X_vol, Z_vol] = vol_element_coords(p)
% returns the coordinates of the center of each volume element in the mesh
X_vol = zeros(p.Zn-1, p.Xn-1);
Z_vol = zeros(p.Zn-1, p.Xn-1);


for j=1:p.Xn-1
    for i =1:p.Zn-1
        X_vol(i,j) = (p.X(i,j) + p.X(i+1,j) + p.X(i,j+1) + p.X(i+1,j+1) ) /4;
        Z_vol(i,j) = (p.Z(i,j) + p.Z(i+1,j) + p.Z(i,j+1) + p.Z(i+1,j+1) ) /4;
    end
end

end

