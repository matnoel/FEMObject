function w = minus(u,v)
% function w = minus(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u)-double(v);
