function w = plus(u,v)
% function w = plus(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u)+double(v);
