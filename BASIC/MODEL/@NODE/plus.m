function u = plus(u,v)

if isa(u,'NODE') && isa(v,'VECTEUR')
    u.POINT = u.POINT + v ; 
elseif isa(v,'NODE') && isa(u,'VECTEUR')
    v.POINT = v.POINT + u ; 
end

    