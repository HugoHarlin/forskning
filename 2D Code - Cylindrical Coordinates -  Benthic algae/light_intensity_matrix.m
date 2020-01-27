function [I] = light_intensity_matrix(p,I,A)
% Calculates the light intensity at all grid points with Beer-Lamberts law,
% using the trapezoidal method to compute the integrals.
% Assumes A is a matrix with the dimensions [p.Yn, p.Xn]
% returns a matrix I with dimensions [p.Yn, p.Xn]

for x = 1:p.Xn-1
    integral = 0;
    for y=1:p.Yn-1
            y_coord = 0.25*(p.Y(y,x) + p.Y(y+1,x) + p.Y(y,x+1) + p.Y(y+1,x+1));
            integral = integral + y_coord *A(y,x);
            I(y,x) = p.I0 * exp(-p.k*integral -p.kbg*y_coord);
    end
end

end

