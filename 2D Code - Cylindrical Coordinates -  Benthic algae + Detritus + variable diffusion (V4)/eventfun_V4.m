function [x,isterm,dir] = eventfun_V4(t,Y,p)
%% Evaluating size of derivatives
t = 0; %the t variable isn't used but needs to be defined.
rhss = rhs_function_V4(t,Y,p);
x = norm(rhss) - 1e-7;
isterm = 1;
dir = -1;

if(0) % seasonal mixing one a year, all
    if(mod(t,365)-0.001 < 0) % mixing once a year
        
        % extracting concentrations
        A  = Y_t(end,1:(p.Xn-1)*(p.Zn-1));
        Rd = Y_t(end,(p.Xn-1)*(p.Zn-1)+1 : 2*(p.Xn-1)*(p.Zn-1));
        D  = Y_t(end,2*(p.Xn-1)*(p.Zn-1)+1 : 3*(p.Xn-1)*(p.Zn-1));
        
        % Reshaping the vectors into matrices
        A  = reshape(A, [p.Xn-1, p.Zn-1]);
        Rd = reshape(Rd, [p.Xn-1, p.Zn-1]);
        D  = reshape(D, [p.Xn-1, p.Zn-1]);
        A  = A';
        Rd = Rd';
        D  = D';
        
        %Calculating nutrient content at t = Tend
        n_algae     = p.volumes_cyl.*A*p.q;
        n_dissolved = p.volumes_cyl.*Rd;
        n_detritus  = p.volumes_cyl.*D;
    end
end

end









