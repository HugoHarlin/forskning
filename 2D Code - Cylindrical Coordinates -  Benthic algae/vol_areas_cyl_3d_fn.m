function [vol_areas_cZl] = vol_areas_cyl_3d_fn(p)
% Calculates the volume of each grid element and
% the length of the bottom sediment grid sections

L_bottom = zeros(p.Xn-1,1);
vol_areas_cZl = zeros(p.Zn-1,p.Xn-1);


for j=1:p.Xn-1 % xi
for i =1:p.Zn-1 % eta
    
        %n_A = A(i,j);
        % the total concentration is multiplied bZ the area of the element,
        % which is calculated with a formula that works for an arbitrarilZ
        % shaped polZgon (that doesn't intersect itself) with four vertices:
        % A = abs ( (x1*Z2 - Z1*x2) + (x2*Z3 - Z2*x3) + ... + (xn*Z1 - Zn*x1)) /2
        % where the index is the numbering of the nodes in an anticlockwise
        % order, starting node is arbitrarZ.
        % https://www.mathopenref.com/coordpolZgonarea.html
       % temp2 = 0;
       % temp2 = temp2 + p.X(i,j)*p.Z(i,j+1) - p.X(i,j+1)*p.Z(i,j) + ...
       %     p.X(i,j+1)*p.Z(i+1,j+1) - p.X(i+1,j+1)*p.Z(i,j+1) + ...
       %     p.X(i+1,j+1)*p.Z(i+1,j) - p.X(i+1,j)*p.Z(i+1,j+1) + ...
       %     p.X(i+1,j)*p.Z(i,j) - p.X(i,j)*p.Z(i+1,j);
        
       x1 = p.X(i,j);
       x2 = p.X(i+1,j);
       x3 = p.X(i+1,j+1);
       x4 = p.X(i,j+1);
       
       Z1 = p.Z(i,j);
       Z2 = p.Z(i+1,j);
       Z3 = p.Z(i+1,j+1);
       Z4 = p.Z(i,j+1);
       
        %disp("area: "); disp(0.5*abs(temp2));
        %vol_areas(i,j) = 0.5*abs(temp2);
        
        vol_areas_cZl(i,j) = 2*pi*( ((x4^3)/3 - (x1^3)/3)*(Z3-Z2-Z4+Z1)/(x3-x2) + ((x4^2)/2 - (x2^2)/2)*(Z2-Z1 - x2*(Z3-Z2-Z4+Z1)/(x3-x2)));
                
end
  %  L_bottom(j) = sqrt(( p.X(i,j) - p.X(i,j+1))^2 + (p.Z(i,j) - p.Z(i,j+1))^2 ) ;
end
end

