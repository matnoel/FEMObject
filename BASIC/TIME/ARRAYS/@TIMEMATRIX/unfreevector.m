function u = unfreevector(u,S,k)
% function u = unfreevector(u,S,k)

if nargin==2
    k=1;
end
u.value = unfreevector(S,u.value,k);
u.s(k)=size(u.value,k);