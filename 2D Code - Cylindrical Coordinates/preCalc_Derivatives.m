function [dX_dXi_preCalc, dX_dEta_preCalc, dZ_dXi_preCalc, dZ_dEta_preCalc] = preCalc_Derivatives(p)
% Evaluation of transform derivatives at nodal values used in numerical
% method. Since the points of evaluation and the derivatives are not time
% dependent this can be done preemtively to save time.


% The third dimension of the matrices correspond to the integral number,
% with the same ordering as the finite volume method comvention used.
dX_dXi_preCalc = 0*ones(p.Zn-1, p.Xn-1,4);
dX_dEta_preCalc = 0*ones(p.Zn-1, p.Xn-1,4);
dZ_dXi_preCalc = 0*ones(p.Zn-1, p.Xn-1,4);
dZ_dEta_preCalc = 0*ones(p.Zn-1, p.Xn-1,4);


% Internal Nodes
for i = 1:p.Zn-1 % eta
    for j = 1:p.Xn-1 % xi
        
        %integral 1
        dX_dXi_preCalc(i,j,1)  = dX_dXi(p,i-0.5,j);
        dX_dEta_preCalc(i,j,1) = dX_dEta(p,i-0.5,j);
        dZ_dXi_preCalc(i,j,1)  = dZ_dXi(p,i-0.5,j);
        dZ_dEta_preCalc(i,j,1) = dZ_dEta(p,i-0.5,j);
        
        dX_dXi_preCalc(i,j,2)  = dX_dXi(p,i-1,j-0.5);
        dX_dEta_preCalc(i,j,2) = dX_dEta(p,i-1,j-0.5);
        dZ_dXi_preCalc(i,j,2)  = dZ_dXi(p,i-1,j-0.5);
        dZ_dEta_preCalc(i,j,2) = dZ_dEta(p,i-1,j-0.5);
        
        dX_dXi_preCalc(i,j,3)  = dX_dXi(p,i-0.5,j-1);
        dX_dEta_preCalc(i,j,3) = dX_dEta(p,i-0.5,j-1);
        dZ_dXi_preCalc(i,j,3)  = dZ_dXi(p,i-0.5,j-1);
        dZ_dEta_preCalc(i,j,3) = dZ_dEta(p,i-0.5,j-1);        
        
        dX_dXi_preCalc(i,j,4)  = dX_dXi(p,i,j-0.5);
        dX_dEta_preCalc(i,j,4) = dX_dEta(p,i,j-0.5);
        dZ_dXi_preCalc(i,j,4)  = dZ_dXi(p,i,j-0.5);
        dZ_dEta_preCalc(i,j,4) = dZ_dEta(p,i,j-0.5);            
        
        % Integral 1 - alternative evaluation of derivatives
%         dAdt(i,j) = dAdt(i,j) +  p.dx * 1/dX_dXi(p,i-0.5,j) * ...
%             ( dY_dEta(p,i-0.5,j) *0.25*(A(i-1,j+1) + 2*A(i,j+1) + A(i+1,j+1) - A(i+1,j) - 2*A(i,j) - A(i-1,j)) ...
%             - dY_dXi(p,i-0.5,j) * 0.25 *( A(i+1,j) + A(i+1,j+1) - A(i-1,j) - A(i-1,j+1) ));
        

        % Integral 2 - alternative derivative approximation
%         dAdt(i,j) = dAdt(i,j) + (p.dx * dY_dXi(p,i-1,j-0.5)/dX_dXi(p,i-1,j-0.5)* ...
%             0.25*(A(i-1,j+1) + A(i,j+1) - A(i-1,j-1) - A(i,j-1)) - 1/( dX_dXi(p,i-1,j-0.5) ...
%         *dY_dEta(p,i-1,j-0.5))*0.25*(A(i,j-1) + 2*A(i,j) + A(i,j+1) - A(i-1,j+1) - 2*A(i-1,j) - A(i-1,j-1) )* ...
%             (p.dx*dY_dXi(p,i-1,j-0.5)^2 + p.dy*dX_dXi(p,i-1,j-0.5)^2));
        

                % Integral 3 - alternative derivative approximation
%         dAdt(i,j) = dAdt(i,j) + sqrt( dX_dEta(p,i-0.5,j-1)^2 + dY_dEta(p,i-0.5,j-1)^2) / abs( dX_dXi(p,i-0.5,j-1)* ...
%             dY_dEta(p,i-0.5,j-1)) * p.dx*( dY_dXi(p,i-0.5,j-1)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i-1,j-1)-A(i-1,j)) - ...
%             dY_dEta(p,i-0.5,j-1) * 0.25*(A(i+1,j) + 2*A(i,j) + A(i-1,j) - A(i-1,j-1) - 2*A(i,j-1) - A(i+1,j-1) ));

        % Integral 4 - alternative derivative approximation
%         dAdt(i,j) = dAdt(i,j) + 1/abs(dX_dXi(p,i,j-0.5) * dY_dEta(p,i,j-0.5))*p.dy*dX_dXi(p,i,j-0.5)^2 * ...
%             0.25*(A(i+1,j-1) + 2*A(i+1,j) + A(i+1,j+1) - A(i,j+1) -2*A(i,j) - A(i,j-1) ) + p.dx*dY_dXi(p,i,j-0.5)^2* ...
%             0.25*(A(i+1,j-1) + 2*A(i+1,j) + A(i+1,j+1) - A(i,j+1) -2*A(i,j) - A(i,j-1) ) -  ...
%             p.dx*dY_dXi(p,i,j-0.5)*dY_dEta(p,i,j-0.5)*0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j-1) - A(i+1,j-1));        

    end
end


end

