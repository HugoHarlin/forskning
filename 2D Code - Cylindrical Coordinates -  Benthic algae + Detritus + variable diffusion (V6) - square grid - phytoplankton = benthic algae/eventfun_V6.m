function [x,isterm,dir] = eventfun_V6(t,Y,p)
%% Evaluating size of derivatives
rhss = rhs_function_V6(t,Y,p);
x = norm(rhss) - 1e-7;

% if(exist('p.seasonal_mixing'))
%     if(p.seasonal_mixing)
%         if(p.seasonal_t > p.mixing_period) % the mixing event happens every p.mixing_period days
%             x = 0;
%         end
%     end
% end

isterm = 1;
dir = -1;

end









