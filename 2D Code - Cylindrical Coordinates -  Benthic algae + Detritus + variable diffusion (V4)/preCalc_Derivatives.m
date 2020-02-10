function [dX_dXi_preCalc, dX_dEta_preCalc, dZ_dXi_preCalc, dZ_dEta_preCalc] = preCalc_Derivatives(p)
% Evaluation of transform derivatives at nodal values used in numerical
% method. Since the points of evaluation and the derivatives are not time
% dependent this can be done preemtivelZ to save time.


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

    end
end


end

