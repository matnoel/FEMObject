function mat = getmaterial(mat,k)
% function mat = getmaterial(mat,k)

if nargin==2 && k~=mat.number
   error('Wrong material number') 
end
