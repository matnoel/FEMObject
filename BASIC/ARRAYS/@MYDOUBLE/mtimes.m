function w = mtimes(u,v)
% function w = mtimes(u,v)

if isa(u,'MYDOUBLE') && isa(v,'MYDOUBLE')
    w = u;
    w.double = u.double*v.double;
elseif isa(u,'MYDOUBLE') 
    w = u;
    w.double = u.double*double(v);
else 
    w = v;
    w.double = double(u)*v.double;
end
