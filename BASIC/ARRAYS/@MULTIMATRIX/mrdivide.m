function w=mrdivide(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'double')
 if isa(u.value,'cell')
     w=u;
     for i=1:numel(u.value)
     w.value{i}=w.value{i}/v;    
     end
 else
     
    if all(size(v)==1)
    w=u;
    w.s=u.s;
    w.value = u.value/v;
    end
 end
 
elseif isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX')
    if isa(u.value,'cell') & isa(v.value,'cell')
    if all(u.sm==v.sm)
    w=u;
     for i=1:numel(u.value)
     w.value{i}=u.value{i}/v.value{i};    
     end
    else
        error('not the same multidimensions') 
    end
    
    else
        error('mrdivide not defined') 
    end
    
else    
       error('mrdivide not defined') 
end

