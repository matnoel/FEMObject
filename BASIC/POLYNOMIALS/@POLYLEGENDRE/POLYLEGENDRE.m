function h = POLYLEGENDRE(varargin)
% function h = POLYLEGENDRE()
% Legendre polynomials Le_n(x)
% polynomes orthonorm�s sur [-1,1] pour la mesure 1/2


h=struct();
param = struct();
domain = [-1,1];
hp = RANDPOLY('legendre',param,domain);
h = class(h,'POLYLEGENDRE',hp);
