function [X,PC] = PCMODEL(a,varargin)
% function X = PCMODEL(a,'order',p)
% a: RANDVARS (elles vont etre decomposees sur le chaos)
% function sto = PCMODEL(a,'order',p,'typebase',typebase)
% typebase=1 : p ordre maxi des polynomes multidimensionnel (par defaut)
% typebase=2 : p ordre maxi des polynomes unidimensionnel
%       p peut etre un vecteur : ordre different par dimension
%
% function sto = PCMODEL(a,'order',p,'randpolys',H)
% H : RANDPOLYS
% H{i} peut etre un polynome classique (POLYLEGENDRE, ...) ou un POLYFE
% pour EF au niveau stochastique, ou un POLYWAVELETS pour des ondelettes
%
% function sto = PCMODEL(a,'order',p,'fedim',i,'femesh',n)
% FE au niveau stochastique pour les dimensions i
% n : nombre de subdivisions de [0,1]
% n : tableau de (length(i)) cellules. n{j} contient le decoupage de [0,1]
% pour la dimension i(j)
%
% function sto = PCMODEL(a,'order',p,'waveletdim',i,'waveletlevel',n)
% Ondelettes au niveau stochastique pour les dimensions i
% n : niveau de resolution (0, 1, ...)
% n : tableau de (length(i)) cellules. n{j} contient le niveau de
% resolution pour la dimension i(j)
%
% function sto = PCMODEL(a,'order',p,'pcgdim',i)
% PC generalise au niveau stochastique pour les dimensions i
% appel de RANDPOLY(a{i}) pour determiner la base polynomiale de la VA a{i}
%
% function sto = PCMODEL(a,'order',p,'pcg')
% PC generalise au niveau stochastique pour toutes les dimensions
% appel de RANDPOLYS(a)
%
% function sto = PCMODEL(a,...,'nomasse')
% don't compute the masse matrix

if length(a)==0
    X = PCMODEL();
else
    disp('---- CREATION OF STOCHASTIC MODEL -----')
    p = getcharin('order',varargin);
    if isempty(p)
        error('preciser l''ordre de la decomposition (degre des polynomes)')
    elseif length(p)==1
        p=repmat(p,1,max(1,a.M));
    end
    tol=getcharin('tol',varargin);
    if ~isempty(tol) && length(tol)==1
        tol=repmat(tol,1,a.M);
    end
    
    % if getcharin('double',varargin)
    %     varargin = setcharin('order',varargin,2*p);
    % end
    
    % H = create_randpolys(a,varargin{:})
    
    varargin1 = setcharin('order',varargin,p);
    H = create_randpolys(a,varargin1{:});
    
    disp('-> Creation of multidimensional chaos')
    typebase=getcharin('typebase',varargin,1);
    
    PC = POLYCHAOS(H,p,'typebase',typebase) ;
    
    if ~ischarin('nomasse',varargin)
        PC = calc_masse(PC);
    end
    
    disp('-> Unidimensional decomposition of random variables')
    X = project(a,PC);
    
    X = PCMODEL(X,varargin{:});
end
