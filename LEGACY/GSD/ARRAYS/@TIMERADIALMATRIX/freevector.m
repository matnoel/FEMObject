function u = freevector(u,S,k)
% function u = freevector(u,S,k)

if nargin==2
    k=1;
end
u.V = freevector(S,u.V,k);
