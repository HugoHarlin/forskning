function [x,isterm,dir] = eventfun(t,Y,p)
   %% Eventfunction stopping the simulation when steady state is reached.
t = 0; %the t variable isn't used but needs to be defined.
rhs = dAdt(t,Y,p);
x = norm(rhs)  - 1e-9;
isterm = 1;
dir = 0;  %or -1, doesn't matter

end







