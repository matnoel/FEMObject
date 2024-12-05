function N = trans(N,v,liste)

if isa(v,'VECTEUR')
    N.POINT=N.POINT+v;
end