function [x,isterm,dir] = eventfun(t,Y,p)
   %% Evaluating size of derivatives
t = 0; %the t variable isn't used but needs to be defined.
rhss = dAdt_efficient_correct(t,Y,p);
x = norm(rhss)  - 1e-8;
isterm = 1;
dir = 0;  %or -1, doesn't matter

end









