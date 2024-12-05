function w = times(u,v)
% function w = times(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u).*double(v);
