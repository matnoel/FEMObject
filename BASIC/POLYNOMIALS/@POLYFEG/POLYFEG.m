function h = POLYFEG(n,rv,varargin)
% function h = POLYFEG(n,rv,varargin)
% base element fini orthonorm�s sur une mesure arbitraire 
% rv : RANDVAR (germe)
%  
% n : si scalaire : decoupage uniforme en n elements
%     si vecteur  : coordonnes des noeuds du maillage
%
% function h = POLYFE(n,r,'node')
% n correspond aux coordonnees des noeuds du maillage
%
% function h = POLYFE(n,r,'elem')
% n : matrice N-by-2 
% chaque ligne correspondant aux noeuds d'un element

if nargin==0
    h=struct();
    h.rv = RANDVAR();
    h = class(h,'POLYFEG',POLYFE());
else    
    if ~isa(rv,'RANDVAR')
        error('rentrer une RANDVAR')
    end
    h.rv = rv;
    h = class(h,'POLYFEG',POLYFE(n,[],varargin{:}));

end
