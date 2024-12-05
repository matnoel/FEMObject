function H = modify(H,varargin)
% function H = modify(H,'fedim',i,'femesh',n)
% FE au niveau stochastique pour les dimensions i 
% n : nombre de subdivisions de [0,1]
% n : tableau de (length(i)) cellules. n{j} contient le decoupage de [0,1]
% pour la dimension i(j)
%
% function H = modify(H,'pcg)
% PC generalise au niveau stochastique 
%
% function H = modify(H,'pcgdim',i)
% PC generalise au niveau stochastique pour les dimensions i
% appel de RANDPOLY(a{i}) pour determiner la base polynomiale de la VA a{i}
%

if ischarin('pcg',varargin)
H=RANDPOLYS(RANDVARS(H)); 
elseif ischarin('pcgdim',varargin)
  pcgdim = getcharin('pcgdim',varargin);  
  H(pcgdim) = RANDPOLYS(RANDVARS(H.h(pcgdim)));
end
fedim = getcharin('fedim',varargin);
if ~isempty(fedim)
n=getcharin('femesh',varargin); 
if ~isa(n,'cell') || length(n)~=length(fedim)
error('femesh n''est pas correct')
end
for k=1:length(fedim)
H(fedim(k)) = POLYFE(n{k},p(k));
end
end
