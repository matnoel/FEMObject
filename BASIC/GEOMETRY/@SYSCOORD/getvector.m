function v = getvector(v,dim)
% function v = getvector(v,dim)

if nargin==1
    v = VECTEUR(v.MYDOUBLEND);
else
    v = VECTEUR(v.MYDOUBLEND(:,dim));
end
