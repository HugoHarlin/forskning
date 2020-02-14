function [x,isterm,dir] = eventfun_V4(t,Y,p)
   %% Evaluating size of derivatives
t = 0; %the t variable isn't used but needs to be defined.
rhss = rhs_function_V4(t,Y,p);
x = norm(rhss) - 1e-7;
isterm = 1;
dir = -1;

end









