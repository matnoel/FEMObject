function V = subsref(V,s)
% function V = subsref(V,s)

if  strcmp(s.type,'()')
    V.dim = V.dim(s.subs{1});
else
    error('ecrire V(i) pour obtenir un sous-PRODUCSPACE')
end






