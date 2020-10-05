function [x,isterm,dir] = eventfun_V6(t,Y,p)
%% Evaluating size of derivatives
rhss = rhs_function_V6(t,Y,p);
x = norm(rhss) - 1e-7;

% if(p.seasonality)
%     if(t(end) > 365) % the mixing event happens every p.mixing_period days
%         if(norm(Y(end,:)) - norm(Y(end-365,:)) < 1e-7)
%             x = 0;
%         end
%     end
% end

isterm = 1;
dir = -1;

end









