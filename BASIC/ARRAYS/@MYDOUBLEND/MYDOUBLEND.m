function u = MYDOUBLEND(v)
% function u = MYDOUBLEND(v)

if nargin==0
    u = MYDOUBLEND([]);
elseif strcmp(class(v),'MYDOUBLEND')
    u = v;
else
    u = struct();
    u.double = full(double(v));
    
    u = class(u,'MYDOUBLEND');
    superiorto('MYDOUBLE')
end
