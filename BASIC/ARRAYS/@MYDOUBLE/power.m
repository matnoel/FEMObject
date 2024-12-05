function w = power(u,v)
% function w = power(u,v)

if isa(u,'MYDOUBLE')
    w = u;
else
    w = v;
end

w.double = double(u).^double(v);
