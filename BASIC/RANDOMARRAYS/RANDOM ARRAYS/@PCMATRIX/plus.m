function w=plus(u,v)

if isa(u,'PCMATRIX') && isa(v,'PCMATRIX')
    w=u;
    w.MULTIMATRIX = u.MULTIMATRIX + v.MULTIMATRIX ;    
elseif isa(u,'PCMATRIX') && isa(v,'PCRADIALMATRIX')
   w = u + PCMATRIX(v) ;
elseif isa(u,'PCRADIALMATRIX') && isa(v,'PCMATRIX')
   w = PCMATRIX(u) + v ;
elseif isa(v,'double')
   w = u + times(v,one(u.POLYCHAOS)); 
elseif isa(u,'double')
   w = v + times(u,one(v.POLYCHAOS));     
else
    class(u)
    class(v)
    
       error('plus not defined') 
end

