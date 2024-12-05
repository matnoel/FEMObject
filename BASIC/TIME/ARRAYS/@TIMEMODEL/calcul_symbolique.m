function M = calcul_symbolique(ordre,p,d)
% function M = calcul_symbolique(ordre,p,d)
% Fiche de calcul symbolique :
%   calcul des integrales int[ u_1*...*u_ordre] sur [0,1]
%   avec u_i choisie parmis les fonction d'interpolation
%   de degre p sur [0,1].
%   Si d est precise, alors calcul de:
%             int[ u_1*...*d(u_d)/dt*...*u_ordre] sur [0,1]

% ordre de la forme :
%    ordre=1   ->  int[ u ]
%    ordre=2   ->  int[ u1*u2 ]
%    ordre=3   ->  int[ u1*u2*u3 ] ...

% degre de l'interpolation
%    p=0   ->  fonction constante (1)
%    p=1   ->  fonction lineaire (2)
%    p=2   ->  fonction quadratique (3)

% derivation pour un certain ordre
%    d=0   ->  aucune derivation
%    d=1   ->  derivation des fonctions de formes de l'ordre 1
%    d=2   ->  derivation des fonctions de formes de l'ordre 2...
if nargin<3
    d=0;
end

if d>ordre
    error('donner un indice d pour la dimenions a deriver');
end
D0=1:ordre;
D1=[];
if d~=0
    D0(d)=[];
    D1=d;
end


% variables symboliques
syms xi PROD;

% construction de la base phi
switch p
    case 0
        phi={1};
    case 1
        phi={1-xi , xi};
    case 2
        phi={2*(xi-1)*(xi-1/2);...
             -4*xi*(xi-1);...
              2*xi*(xi-1/2) };
    otherwise
        error('pas implemente...')
end

% calcul des integrales :
s=length(phi)*ones(1,ordre);
M=zeros(s);
if ordre==1, M=zeros(s,1);end
sub=cell(1,ordre);
for i=1:numel(M)
    [sub{:}]=ind2sub(s,i);
    PROD=1;
    for o=D0
        PROD=PROD*phi{sub{o}};
    end
    for o=D1
        PROD=PROD*diff(phi{sub{o}});
    end
    M(i)=double(int(PROD,0,1));
end


% calcul evaluation sur le bord de [0,1]
%...











