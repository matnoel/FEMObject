function w = mtimes(u,v)
% function w = mtimes(u,v)

if isa(u,'MYDOUBLEND')
    w = u;
else
    w = v;
end

[u,v] = samesizeND(u,v);
w.double = times3D(reshape3D(u),reshape3D(v));
w.double = reshape(w.double,[size2D(w),sizeND(u)]);
