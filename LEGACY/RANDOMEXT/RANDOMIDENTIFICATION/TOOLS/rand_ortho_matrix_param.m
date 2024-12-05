function phi = rand_ortho_matrix_param(n,p,lb,ub)
% phi = rand_ortho_matrix_param(n,p)
% tirage aleatoire des parametres d'une matrice orthogonale n-by-p
% phi = rand_ortho_matrix_param(n,p)

dimvariete = n*p-n*(n+1)/2;
if nargin==2
    [lb,ub]=ortho_matrix_param_bounds(n,p);
else
    if numel(lb)==1
        lb = repmat(lb,dimvariete,1);
    end   
    if numel(ub)==1
        ub = repmat(ub,dimvariete,1);
    end   
end

phi = lb+(ub-lb).*rand(dimvariete,1);


