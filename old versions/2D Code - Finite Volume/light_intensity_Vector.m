function [I] = light_intensity_Vector(p,I,A)
% Calculates the light intensity at all grid points with Beer-Lamberts law,
% using the trapezoidal method to compute the integrals.
% Assumes A is a matrix with the dimensions [p.Yn, p.Xn]
% returns a matrix I with dimensions [p.Yn, p.Xn]

for x = 1:p.Xn
    % disp("x: "); x
   
    for y=1:p.Yn
        % disp("y: "); y
        integral = 0;
        if(y == 1)
            I(y,x) = p.I0;
        else
            for i=2:y
                integral = integral + (p.Y(i,x)+p.Y(i-1,x))/2 *(A(i,x)+A(i-1,x))/2;
            end
            I(y,x) = p.I0 * exp(-p.k*integral -p.kbg*p.Y(y,x));
        end
    end
end

end

