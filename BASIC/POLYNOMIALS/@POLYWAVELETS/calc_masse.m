function masse=calc_masse(h,p,p2,varargin)
% masse=calc_masse(h,p,p2)
%
% calcul des matrice de masse 
% masse : MULTIMATRIX (p2+1)-by-(p2+1) et de longueur (p+1) 
% E(hk * a * b) = a'*masse{k+1}*b  ou a et b sont les vecteurs
% representant les coeff de a et b sur la base (h0, ... , hp2) 

if nargin==2
p2=p;
elseif p2~=2
    error('non prevu')
end
if ischarin('display',varargin)
fprintf(['Calcul des masse pour ' get(h,'type') ' ... '])
end

if p~=getparam(h,'p')
    error('degree of wavelet does not coincide')
end

error('a programmer')

masse=MULTIMATRIX(masse,[p2+1,p2+1],[p+1,1]);

if ischarin('display',varargin)
    fprintf('\n')
end

if ischarin('double',varargin) || ischarin('cell2mat',varargin)
   masse = cell2mat(masse); 
end
if ischarin('double',varargin)
    warning('remplacer argument double par cell2mat')
end