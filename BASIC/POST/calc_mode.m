function varargout = calc_mode(A,B,varargin)
% function [V,D] = calc_mode(A,B,m)
% calcul les valeurs propres l et vecteurs propres U generalises du systeme
%  A V = l Bs V
% ou Bs est la partie symmetrique de B
% 
% m designe la liste des modes a calculer (numerotes de la plus grande 
%    valeur propre a la plus petite)
% V est matrice dont les colonnes contient les vecteurs propres m
% D est une matrice diagonale dont les termes diagonaux sont les valeurs
% propres m
% 
% function D = calc_mode(A,B,m)
% calcule uniquement les valeurs propres listees dans m

B = (B+B')/2;
m=varargin{1};

if nargout ==1
    D = eigs(A,B,max(m),'SM');
    D = sort(D);
    varargout{1} = D(m);
elseif nargout==2
    [V,D]= eigs(A,B,max(m),'SM');
    [D,renum] = sort(diag(D));
    V = V(:,renum);
    V = V(:,m);
    D = diag(D(m));
    varargout{1} = V;
    varargout{2} = D;
end
