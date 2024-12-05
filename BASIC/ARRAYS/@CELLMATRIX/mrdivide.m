function w=mrdivide(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'double')
 
    if all(size(v)==1)
    w=u;
    w.s=u.s;
    w.value = u.value/v;
    end
    
else
       error('mrdivide not defined') 
end

