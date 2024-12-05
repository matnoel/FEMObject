function w = ldivide(u,v)
% function w = ldivide(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u).\double(v);
