function [dZ_dEta] = dZ_dEta(p,eta,xi)
dZ_dEta = 1/(p.Zn-1) *(p.Lmax + (p.Lmin - p.Lmax) * ( xi/(p.Xn-1)) ^p.alpha);
% dY_dEta = 1/(p.Yn) *(p.Lmax + (p.Lmin - p.Lmax) * ( xi/(p.Xn)) ^p.alpha);
end

