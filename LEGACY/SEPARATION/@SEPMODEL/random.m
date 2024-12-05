function [uD,tirage,newSM] = random(SM,u,dim)
% function [uD,tirage,newSM] = random(SM,u,dim)
% Tirage aleatoire de u :
% 1 - Generer l'operateur d'evaluation 'tirage' (Id sur les dimensions
%     deterministes et un vecteur rand()' sur les autres)
% 2 - tirage*u -> deterministe inchangees
%              -> produit scalaire sur stochastique
% 3 - Eliminer les dimensions devenues scalaires

if nargin==2
    dim = 1:SM.dim;
end

% Generation du tirage
tirage = cell(1,SM.dim);
s = size(u);
for d=1:SM.dim
    if strcmp(SM.F{d}.type,'STOCH')&& ~isempty(find(dim==d, 1))
        % On est sur une dimension stoch
        tirage{d} = full(random(SM.F{d}.model));
    else
        % Ne pas toucher a la dimension :
        tirage{d} = speye(s(d));
    end
end
tirage = SEPMATRIX(tirage);
if isa(u,'HSEP')
    tirage = HSEPMATRIX(tirage,u.tree);
end

% Finaliser :
if min(size(u.F{1,1}))==1
    % Cas ou u est un vecteur
    uD = scalexp(tirage*u);
else
    % Cas ou u est un operateur
    uD = scalexp((tirage*u)*tirage');
end

% arguments de sortie supplementaires :
if nargout>1
    S(1).type = '()';
    S(1).subs = {spatialdim(SM)};
    newSM = subsref(SM,S);
end
