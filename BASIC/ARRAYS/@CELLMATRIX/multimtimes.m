function w=multimtimes(u,v)
%warning('nouvelle fonction')
if isa(u,'MULTIMATRIX')
    u = switchmulti(u);
end
if isa(v,'MULTIMATRIX')
    v = switchmulti(v);
end

w=switchmulti(mtimes(u,v));
%if isa(u,'MULTIMATRIX') & isa(v,'double')
%    w=u;
%    w.value = u.value * v ;
%    w.sm = [size(w.value,2),1];
%end
