function w = mtimes(u,v)
% function w = mtimes(u,v)

if isa(u,'FENODEFIELD') && isa(v,'FENODEFIELD')
    w = u;
    w.value = mtimes(u.value,v.value);
elseif isa(u,'FENODEFIELD') 
    w = u;
    w.value = mtimes(u.value,v);
elseif isa(v,'FENODEFIELD') 
    w = v;
    w.value = mtimes(u,v.value);
end

w.size = size(w.value);
