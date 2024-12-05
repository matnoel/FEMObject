function [P,numnode1,numnode2] = calc_P_transfer(M1,M2,ddls)
% function [P,numnode1,numnode2] = calc_P_transfer(M1,M2,ddls)
% M1, M2 : MODEL
% ddls : liste des ddls a bloquer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL 
% ddls = 'all' par defaut
% calcul de la matrice permettant de passer des ddl de M1 aux ddl de M2
% si v1 est un vecteur de ddl defini sur M1
% v2 = P*v1 est le vecteur de ddl correspondant defini sur M2

if nargin<3
    [P,numnode1,numnode2] = calc_P(M1,M2);
else
    [P,numnode1,numnode2] = calc_P(M1,M2,ddls);
end
