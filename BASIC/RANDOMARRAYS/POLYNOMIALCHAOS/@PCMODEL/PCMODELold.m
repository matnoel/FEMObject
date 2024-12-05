function sto = PCMODEL(a,varargin)
% function X = PCMODEL(a,'order',p)
% a: RANDVARS (elles vont etre decomposees sur le chaos)
%  function sto = PCMODEL(a,'order',p,'typebase',typebase)
%  typebase=1 : p ordre maxi des polynomes multidimensionnel (par defaut)
%  typebase=2 : p ordre maxi des polynomes unidimensionnel
%       p peut etre un vecteur : ordre different par dimension
%
% function sto = PCMODEL(a,'order',p,'randpolys',H)
% H : RANDPOLYS
% H{i} peut etre un polynome classique ou un POLYFE pour EF au niveau stochastique
%
% function sto = PCMODEL(a,'order',p,'fedim',i,'femesh',n)
% FE au niveau stochastique pour les dimensions i
% n : nombre de subdivisions de [0,1]
% n : tableau de (length(i)) cellules. n{j} contient le decoupage de [0,1]
% pour la dimension i(j)
%
% function sto = PCMODEL(a,'order',p,'pcgdim',i)
% PC generalise au niveau stochastique pour les dimensions i
% appel de RANDPOLY(a{i}) pour determiner la base polynomiale de la VA a{i}
%
% function sto = PCMODEL(a,'order',p,'pcg')
% PC generalise au niveau stochastique pour toutes les dimensions
% appel de RANDPOLYS(a)
%
% function sto = PCMODEL(a,'double')
% decomposition des variables sur des polynomes de degre 2*p
% calcul de la matrice masse en consequence

if nargin==1  && isa(a,'PCMODEL')
    sto=a;
elseif nargin==1 && isa(a,'PCMATRIX')
    PC = getPC(a);
    sto.X = a(:);
    sto=class(sto,'PCMODEL',PC);
    superiorto('POLYCHAOS')
else
    error('mauvais arguments')
end

return




a=RANDVARS(a);
disp('---- CREATION OF STOCHASTIC MODEL -----')

H = getclassin('RANDPOLYS',varargin);

if ~isempty(H)
    if length(H)~=a.M
        error('RANDPOLYS doit contenir autant de polynomes que le nombre de VA')
    end
else
    H=RANDPOLYS();
    H(1:a.M)=POLYHERMITE();
    if ischarin('pcg',varargin)
        H=RANDPOLYS(a);
    elseif ischarin('pcgdim',varargin)
        pcgdim = getcharin('pcgdim',varargin);
        H(pcgdim) = RANDPOLYS(a.RV(pcgdim));
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
end

H=setnumber(H,getnumber(a));


p=getcharin('order',varargin);
if isempty(p)
    error('preciser l''ordre de la decomposition (degre des polynomes)')
end

if getcharin('double',varargin)
    varargin = setcharin('order',varargin,2*p);
end

disp('-> Unidimensional decomposition of random variables')
[X,h] = decomppc(a,varargin{:});

disp('-> Creation of multidimensional chaos')
typebase=getcharin('typebase',varargin,1);

if ischarin('double',varargin)
    PC = POLYCHAOS(h,2*p,'typebase',typebase);
    PC2 = POLYCHAOS(h,p,'typebase',typebase);
    PC = calc_masse(PC,PC2);
else
    PC = POLYCHAOS(h,p,'typebase',typebase);
    PC = calc_masse(PC);
end

sto =struct();

for i=1:a.M
    sto.X{i}=project(X{i},PC);
end
disp(' ')


sto.X = vertcat(sto.X{:});
sto.RANDVARS = a;
sto=class(sto,'PCMODEL',PC);
superiorto('POLYCHAOS')
