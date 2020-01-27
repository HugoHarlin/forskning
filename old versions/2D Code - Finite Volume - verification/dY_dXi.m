function [dY_dXi] = dY_dXi(p,eta,xi)
dY_dXi =  p.alpha*eta/((p.Yn-1)*(p.Xn-1)) * (p.Lmin - p.Lmax) * (xi/(p.Xn-1)) ^(p.alpha-1) ;
% dY_dXi =  p.alpha*eta/((p.Yn)*(p.Xn)) * (p.Lmin - p.Lmax) * (xi/(p.Xn)) ^(p.alpha-1) ;
end

