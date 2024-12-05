function u = freevector(u,S,k)
% function u = freevector(u,S,k)

if nargin==2
    k=1;
end
u.value = freevector(S,u.value,k);
u.s(k)=size(u.value,k);