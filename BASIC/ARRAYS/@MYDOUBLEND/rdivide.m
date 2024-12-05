function w = rdivide(u,v)
% function w = rdivide(u,v)

if isa(u,'MYDOUBLEND')
    w = u;
else
    w = v;
end

[u,v] = samesize(u,v);
w.double = u./v;
