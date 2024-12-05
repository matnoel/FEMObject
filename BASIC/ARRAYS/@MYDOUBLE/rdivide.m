function w = rdivide(u,v)
% function w = rdivide(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u)./double(v);
