function w = plus(u,v)
% function w = plus(u,v)

if isa(u,'MYDOUBLEND')
    w = u;
else
    w = v;
end

[u,v] = samesize(u,v);
w.double = u+v;
