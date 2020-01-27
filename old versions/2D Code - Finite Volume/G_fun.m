function [G] = G_fun(p,A,Rd)
% Calculates teh algal production G at each point in the computation grid
% and returns it as a matix

% Calculation of the light intensity at each grid element center
I = zeros(p.Yn-1,p.Xn-1);
for x = 1:p.Xn-1
    integral = 0;
    for y=1:p.Yn-1
        %y_coord = 0.25*(p.Y(y,x) + p.Y(y+1,x) + p.Y(y,x+1) + p.Y(y+1,x+1));
        integral = integral + p.Y_vol(y,x) *A(y,x);
        I(y,x) = p.I0 * exp(-p.k*integral -p.kbg*p.Y_vol(y,x));
    end
end

% Calculation of the algal production G(R,I)
G =  p.Gmax .* Rd./(Rd+p.M) .* I./(I+p.H);
end

