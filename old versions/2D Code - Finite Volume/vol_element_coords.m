function [X_vol, Y_vol] = vol_element_coords(p)
% returns the coordinates of the center of each volume element in the mesh
X_vol = zeros(p.Yn-1, p.Xn-1);
Y_vol = zeros(p.Yn-1, p.Xn-1);


for j=1:p.Xn-1
    for i =1:p.Yn-1
        X_vol(i,j) = (p.X(i,j) + p.X(i+1,j) + p.X(i,j+1) + p.X(i+1,j+1) ) /4;
        Y_vol(i,j) = (p.Y(i,j) + p.Y(i+1,j) + p.Y(i,j+1) + p.Y(i+1,j+1) ) /4;
    end
end

end

