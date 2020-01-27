function [I] = light_intensity(p,A,x,y)
% Calculates the light intensity at the point (x,y) with Beer-Lamberts law,
% using the trapezoidal method to compute the integral.

integral = 0;
if(y == 1)
    I = p.I0;
else
    for i=2:y
        integral = integral + (p.Y(i,x)+p.Y(i-1,x))/2 *(A(i,x)+A(i-1,x))/2;
    end
    I = p.I0 * exp(-p.k*integral -p.kbg*p.Y(y,x));
end

end

