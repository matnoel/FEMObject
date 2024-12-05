function [lb,ub] = ortho_matrix_param_bounds(n,p,varargin)
% borne des parametres d'une variete de Stiefel 
% compacte (ensemble des matrices orthogonales de dimension n-by-p)
% %  cas n=1 : vecteur orthonorme dans R^p
%             les parametres sont les parametres angulaires d'une
%             hypersphere
%             phi = (t1,t2,...,tp-1) in [0,pi]x[0,pi]x...x[-pi,pi]
%
global generalstiefel
if isempty(generalstiefel)
    generalstiefel=0;
end

if n==1 && ~generalstiefel
    lb=zeros(p-1,1);ub=zeros(p-1,1);
    if ~ischarin('symmetry',varargin)  
        ub(1:end)=pi;
        lb(end)=-pi;
    else
        ub(:)=pi;
    end

else
    d=n*p-n*(n+1)/2;
    lb=zeros(d,1);ub=zeros(d,1);
    lb(:)=-pi;ub(:)=pi;
% si angle < 0 : matrice orthogonale negative
% si angle > 0 : matrice orthogonale positive
end

