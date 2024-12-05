function w = mpower(u,v)
% function w = mpower(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u)^double(v);
