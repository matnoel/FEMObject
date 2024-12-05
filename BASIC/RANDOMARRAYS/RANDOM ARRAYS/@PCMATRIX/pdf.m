function [p,x] = pdf(u,x,N,varargin)
% function p = pdf(u,x,N,'ksdensity')
% evalue la pdf de u (PCMATRIX) en x (double)
% N : nb de tirages pour l'evaluation de la pdf
% 'ksdensity' (optionnel)

if nargin<=2 || isempty(N)
N=1e5;
end
us=randomlimit(u,N);

withks = nargin>=4 && (ischarin('ksdensity',varargin) || ischarin('ks',varargin));
if withks
[p,x] = ksdensity(us,x);
else
[p,x] = pdfsample(us,x);
end

