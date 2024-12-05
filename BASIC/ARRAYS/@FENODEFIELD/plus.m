function w = plus(u,v)
% function w = plus(u,v)

if isa(u,'FENODEFIELD') && isa(v,'FENODEFIELD')
    w = u;
    w.value = u.value+v.value;
elseif isa(u,'FENODEFIELD') 
    w = u;
    w.value = u.value+v;
elseif isa(v,'FENODEFIELD') 
    w = v;
    w.value = u+v.value;
end

w.size = size(w.value);
