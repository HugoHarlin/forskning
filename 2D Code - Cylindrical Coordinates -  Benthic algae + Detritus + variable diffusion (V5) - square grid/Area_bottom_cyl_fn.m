function [L_bottom_cyl] = Area_bottom_cyl_fn(p)
% Calculates areas of bottom segments given cylindrical coordiantes.
% Each bottom segment is a horizontal slice of a cone.
L_bottom_cyl = zeros(1,p.Xn-1);
bottom_index = p.Zn;
for i = 1:p.Xn-1
    z1 = p.Z(bottom_index-1,i+1);
    z2 = p.Z(bottom_index,i); 
    
    x1 = p.X(bottom_index-1,i+1);
    x2 = p.X(bottom_index,i);
    
    L_bottom_cyl(i) = 2*pi*sqrt(1+((z1-z2)/(x1-x2))^2)*0.5*(x1^2-x2^2);
    bottom_index = bottom_index -1;
end
end
