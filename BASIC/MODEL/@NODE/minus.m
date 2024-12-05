function u = minus(u,v)

if isa(u,'NODE') & isa(v,'VECTEUR')
    u.POINT = u.POINT - v ; 
elseif isa(v,'NODE') & isa(u,'VECTEUR')
    v.POINT = u  -  v.POINT  ; 
end

    