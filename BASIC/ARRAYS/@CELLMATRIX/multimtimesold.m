function w=multimtimesold(u,v)
warning('utiliser la nouvelle fonction multimtimes')
if isa(u,'MULTIMATRIX') & isa(v,'double')
    w=u;
    w.value = u.value * v ;
    w.sm = [size(w.value,2),1];
end
