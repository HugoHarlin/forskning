function [dZ_dXi] = dZ_dXi(p,eta,xi)
dZ_dXi =  p.alpha*eta/((p.Zn-1)*(p.Xn-1)) * (p.Lmin - p.Lmax) * (xi/(p.Xn-1)) ^(p.alpha-1) ;
% dZ_dXi =  p.alpha*eta/((p.Zn)*(p.Xn)) * (p.Lmin - p.Lmax) * (xi/(p.Xn)) ^(p.alpha-1) ;
end

