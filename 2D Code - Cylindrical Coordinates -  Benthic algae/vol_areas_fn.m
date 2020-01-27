function [vol_areas,L_bottom] = vol_areas_fn(p)
% Calculates the volume of each grid element and
% the length of the bottom sediment grid sections

L_bottom = zeros(p.Xn-1,1);
vol_areas = zeros(p.Yn-1,p.Xn-1);


for j=1:p.Xn-1 % xi
for i =1:p.Yn-1 % eta
    
        %n_A = A(i,j);
        % the total concentration is multiplied by the area of the element,
        % which is calculated with a formula that works for an arbitrarily
        % shaped polygon (that doesn't intersect itself) with four vertices:
        % A = abs ( (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + ... + (xn*y1 - yn*x1)) /2
        % where the index is the numbering of the nodes in an anticlockwise
        % order, starting node is arbitrary.
        % https://www.mathopenref.com/coordpolygonarea.html
        temp2 = 0;
        temp2 = temp2 + p.X(i,j)*p.Y(i,j+1) - p.X(i,j+1)*p.Y(i,j) + ...
            p.X(i,j+1)*p.Y(i+1,j+1) - p.X(i+1,j+1)*p.Y(i,j+1) + ...
            p.X(i+1,j+1)*p.Y(i+1,j) - p.X(i+1,j)*p.Y(i+1,j+1) + ...
            p.X(i+1,j)*p.Y(i,j) - p.X(i,j)*p.Y(i+1,j);
        
        %disp("area: "); disp(0.5*abs(temp2));
        vol_areas(i,j) = 0.5*abs(temp2);
        
end
    L_bottom(j) = sqrt(( p.X(i,j) - p.X(i,j+1))^2 + (p.Y(i,j) - p.Y(i,j+1))^2 ) ;
end
end

