function w = times(u,v)
% function w = times(u,v)

if isa(u,'FENODEFIELD') && isa(v,'FENODEFIELD')
    w = u;
    w.value = times(u.value,v.value);
elseif isa(u,'FENODEFIELD') 
    w = u;
    w.value = times(u.value,v);
elseif isa(v,'FENODEFIELD') 
    w = v;
    w.value = times(u,v.value);
end

w.size = size(w.value);
