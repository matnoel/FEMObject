function n = norm(u,varargin)
% function n = norm(u,varargin)

% Pour deux vecteurs pour l'instant
% n = sqrt(expand(mtimes(u',u)));

% Pour deux matrices (ou vecteurs simples)
n = sqrt(abs(fastprodscal(u,u)));
% n = sqrt(abs(expand(u'*u)));
