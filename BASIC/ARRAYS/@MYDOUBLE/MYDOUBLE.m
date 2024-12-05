function u = MYDOUBLE(v)
% function u = MYDOUBLE(v)

if nargin==0
    v = [];
end
u.double = double(v);
u = class(u,'MYDOUBLE');
