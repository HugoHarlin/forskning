function [vol_areas_cyl] = vol_areas_cyl_3d_fn(p)
% Calculates the volume of each grid element and
% the length of the bottom sediment grid sections

L_bottom = zeros(p.Xn-1,1);
vol_areas_cyl = zeros(p.Zn-1,p.Xn-1);


for j=1:p.Xn-1 % xi
for i =1:p.Zn-1 % eta
    
        %n_A = A(i,j);
        % the total concentration is multiplied by the area of the element,
        % which is calculated with a formula that works for an arbitrarily
        % shaped polygon (that doesn't intersect itself) with four vertices:
        % A = abs ( (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + ... + (xn*y1 - yn*x1)) /2
        % where the index is the numbering of the nodes in an anticlockwise
        % order, starting node is arbitrary.
        % https://www.mathopenref.com/coordpolygonarea.html
       % temp2 = 0;
       % temp2 = temp2 + p.X(i,j)*p.Y(i,j+1) - p.X(i,j+1)*p.Y(i,j) + ...
       %     p.X(i,j+1)*p.Y(i+1,j+1) - p.X(i+1,j+1)*p.Y(i,j+1) + ...
       %     p.X(i+1,j+1)*p.Y(i+1,j) - p.X(i+1,j)*p.Y(i+1,j+1) + ...
       %     p.X(i+1,j)*p.Y(i,j) - p.X(i,j)*p.Y(i+1,j);
        
       x1 = p.X(i,j);
       x2 = p.X(i+1,j);
       x3 = p.X(i+1,j+1);
       x4 = p.X(i,j+1);
       
       y1 = p.Z(i,j);
       y2 = p.Z(i+1,j);
       y3 = p.Z(i+1,j+1);
       y4 = p.Z(i,j+1);
       
        %disp("area: "); disp(0.5*abs(temp2));
        %vol_areas(i,j) = 0.5*abs(temp2);
        
        vol_areas_cyl(i,j) = 2*pi*( ((x4^3)/3 - (x1^3)/3)*(y3-y2-y4+y1)/(x3-x2) + ((x4^2)/2 - (x2^2)/2)*(y2-y1 - x2*(y3-y2-y4+y1)/(x3-x2)));
                
end
  %  L_bottom(j) = sqrt(( p.X(i,j) - p.X(i,j+1))^2 + (p.Y(i,j) - p.Y(i,j+1))^2 ) ;
end
end

