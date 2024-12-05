function u = unfreevector(S,u,k)
% function u = unfreevector(S,u,k)

if nargin==2
    k=1;
end
u.V = unfreevector(S,u.V,k);
