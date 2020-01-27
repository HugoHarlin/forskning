function [dX_dXi] = dX_dXi(p,eta,xi)
 dX_dXi = p.W/(p.Xn-1);
% dX_dXi = p.W/(p.Xn);
end