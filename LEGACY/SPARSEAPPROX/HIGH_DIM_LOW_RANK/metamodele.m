function Y=metamodele(X)

% X données d'entrée en dimension 4 comprises entre -1 et 1 (lois
% uniformes)

% Y sortie correspondante à X


load multi_indices.mat

if size(X,2)~=4
    error('Le nombre de colonne de X doit être égal à 4')
end

sup=isempty(find(X>1,1));
inf=isempty(find(X<-1,1));

if (sup==0 || inf==0)
    error('Les données d''entrée doivent être comprises entre -1 et 1')
end

PSI=Legendre(X,alphaf);

Y=PSI*betaf; 

end


function [ Psi ] = Legendre( X,alpha )

%On g�n�re ici la matrice PSI contenant les polyn�mes de Legendre
%
%entr�es :
% X = matrice contenant en colonne les variables d'entr�e de ton probl�me
% (ici ce sont n�cessairement des U([-1;1]) ind�pendantes!) et en lignes
% les observations
% alpha = la matrice des multi-degr�s (M lignes, nombre de colonnes
% d�pend de la norme choisie) conserv�s apr�s troncature
%
%sorties :
% Psi = matrice de taille (nombre_observation)*(dimension de la base)
% Psi(i,j)=�valuation du j�me polynome (multi-dimensionnel) de la base
% pour la i�me observation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=size(alpha,1);
eval=Legendre_unidim(X,max(alpha(:)));
eval=permute(eval,[1 3 2]);

for i=1:M
    Mat(:,:,i)=eval(:,alpha(i,:)+1,i);
end

Psi=prod(Mat,3);
end

function eval=Legendre_unidim(X,degre_max)
value=X;
x=zeros(size(X,1),size(X,2),degre_max+1);


x(:,:,1)= ones(size(X));
x(:,:,2)= value;

for i=2:degre_max
    % relation de r�currence utilis�e pour obtenir les polyn�mes de
    % Legendre.

   x(:,:,i+1)=(1/i)*((2*(i-1)+1).*value.*x(:,:,i)-(i-1).*x(:,:,i-1));
   
end
eval=x; 

end