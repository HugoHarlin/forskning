function [n_algae, n_dissolved, n_sediment, vol_areas, L_bottom] = Nutrient_content(A,Rd,Rs,p)

%% Total nutrient density of the system at t=0 (to be conserved)

n_sediment = zeros(1,p.Xn-1);
n_algae = zeros(p.Yn-1,p.Xn-1);
n_dissolved = zeros(p.Yn-1,p.Xn-1);
vol_areas = zeros(p.Yn-1,p.Xn-1);
L_bottom = zeros(p.Xn-1,1);
% calculating the totalt nutrient content in the water mass,
% using the trapezoidal method. The average of the four nodes of each
% element in the mesh is calculated, and multiplied by the area of the
% element.
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
        n_A = A(i,j) *p.q* vol_areas(i,j);
        n_Rd = Rd(i,j) * vol_areas(i,j);
        
        n_algae(i,j) = n_A;
        n_dissolved(i,j) = n_Rd;
        n_dissolved(i,j) = n_Rd;
end
    L_bottom(j) = sqrt(( p.X(i,j) - p.X(i,j+1))^2 + (p.Y(i,j) - p.Y(i,j+1))^2 ) ;
end
n_sediment(:) = Rs'.*L_bottom; 

end

