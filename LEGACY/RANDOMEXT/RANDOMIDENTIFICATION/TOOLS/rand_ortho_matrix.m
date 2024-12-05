function A = rand_ortho_matrix(n,p,varargin)
% A = rand_ortho_matrix(n,p)
% tirage aleatoire d'une matrice orthogonale n-by-p
% definition non classique de matrice orthogonale  A*A' = I
%
%  matrice orthogonale parametree (parametrage de la variete de Stiefel)
%  tirage aleatoire des parametres et evaluation de A
%
% A = rand_ortho_matrix(n,p,B)
% tirage aleatoire de A avec parametrisation locale autour de 
% la matrice orthogonale B
%                A = A*B;
%
% A = rand_ortho_matrix(n,p,[],'noparam')
%  tirage aleatoire sans parametrage
%         ->  tirage aleatoire d'une matrice A0 avec composantes dans (-1,1)
%         ->  factorisation de cholesky A0*A0'=R'*R
%         ->  A = R'\A0;
%
% A = rand_ortho_matrix(n,p,B,'noparam')
% tirage local non parametree autour de B

if ischarin('noparam',varargin)
    A = 2*rand(n,p)-1;
    R = chol(A*A');
    A = (R')\A;
    if nargin==4
        A = A*B;    
    end

else
    phi = rand_ortho_matrix_param(n,p);
    A = eval_ortho_matrix_param(n,p,phi,varargin{:});
end
return
