function s = sobol_indices_group(u,g)
% function s = sobol_indices_group(u,g)
% compute the first order (partial) sensitivity indices (i.e. the closed sensivity index)
% with respect to a group of variables g

% PC = getPC(u);
var = expecttimes(u,u)-expect(u).^2;
mo = expect(u);

mi = expectnodim(g,u);
s = (expecttimes(mi,mi)-mo.^2)./var;
