function mat = randomeval(mat,varargin)
% function mat = randomeval(mat,x,RV)
% calcul de realisations des MATERIALS
% on applique randomeval(.,x,RV) a tous les MATERIAL
% 
%  See also MATERIAL/randomeval

for k=1:mat.n
    mat.MAT{k} = randomeval(mat.MAT{k},varargin{:});
end

