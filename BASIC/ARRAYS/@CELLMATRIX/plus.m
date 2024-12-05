function w=plus(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX')
    w=u;
    w.value = u.value + v.value ; 
    
elseif isa(u,'MULTIMATRIX') & isa(v,'double')
    
    if all(size(v)==1)
    value = u.value + v ;
    else
    v=repmat(v(:),1,size(u.value,2));
    value = u.value + v ;
    end
    
    w=u;
    w.s=u.s;
    w.value = value;
    
elseif isa(v,'MULTIMATRIX') & isa(u,'double')

    if all(size(u)==1)
    value = v.value + u ;
    else
    u=repmat(u(:),1,size(u.value,2));
    value = u.value + v ;
    end
    
    w=v;
    w.s=v.s;
    w.value = value;
    
else
       error('plus not defined') 
end

