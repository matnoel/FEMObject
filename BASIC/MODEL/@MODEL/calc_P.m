function [P,numnode1,numnode2] = calc_P(M1,M2,ddls)
% function [P,numnode1,numnode2] = calc_P(M1,M2,ddls)
% M1, M2 : MODEL
% ddls : liste des ddls a bloquer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL 
% ddls = 'all' par defaut
% calcul de la matrice permettant de passer des ddl de M1 aux ddl de M2
% si v1 est un vecteur de ddl defini sur M1
% v2 = P*v1 est le vecteur de ddl correspondant defini sur M2

if nargin<3 || isempty(ddls) || strcmp(ddls,'all')
    ddls = 'all';
elseif isa(ddls,'char')
    ddls = {ddls};
end

P = sparse(M2.nbddl,M1.nbddl);

[~,numnode1,numnode2] = intersect(M1.node,M2.node);
ddl1 = findddl(M1,ddls,numnode1);
ddl2 = findddl(M2,ddls,numnode2);
% ind = sub2ind(size(P),ddl2,ddl1);
% P(ind) = 1;
n1 = length(ddl1);
n2 = length(ddl2);
P(ddl2,ddl1) = speye(n1,n2);
