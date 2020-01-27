function [dY_dEta] = dY_dEta(p,eta,xi)
dY_dEta = 1/(p.Yn-1) *(p.Lmax + (p.Lmin - p.Lmax) * ( xi/(p.Xn-1)) ^p.alpha);
% dY_dEta = 1/(p.Yn) *(p.Lmax + (p.Lmin - p.Lmax) * ( xi/(p.Xn)) ^p.alpha);
end

