function T = fulltsep(v)

if isa(v,'double')
    v=tensor(v);
end


szv=size(v);

F=cell(ndims(v),1);
for d=1:ndims(v)
    F{d}=speye(szv(d));
end

T=splitvectors(TSEPMATRIX(F,v));

