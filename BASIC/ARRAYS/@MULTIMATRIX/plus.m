function w=plus(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX')
    w=u;
    if isa(u.value,'double') & isa(v.value,'double')
    w.value = u.value + v.value ; 
    elseif isa(u.value,'cell') & isa(v.value,'cell')
    if all(u.sm==v.sm)
        w=u;
        for k=1:numel(u.value)
        w.value{k}=u.value{k}+v.value{k};    
        end
        w.s = max(u.s,v.s);
    else
       error('not the same multidimensions') 
    end
    elseif isa(u.value,'double') 
     w = plus(mat2cell(u),v);   
    elseif isa(v.value,'double') 
     w = plus(u,mat2cell(v));             
    end
        
elseif isa(u,'MULTIMATRIX') & isa(v,'double')
    if isa(u.value,'cell')
    w=u;
     for k=1:numel(u.value)
        w.value{k}=u.value{k}+v;    
     end
     w.s = max(u.s,size(v));   
    else
    if all(size(v)==1)
    value = u.value + v ;
    else
    v=repmat(v(:),1,size(u.value,2));
    value = u.value + v ;
    end
    
    w=u;
    w.s=u.s;
    w.value = value;
    end    
elseif isa(v,'MULTIMATRIX') & isa(u,'double')
    if isa(v.value,'cell')
    w=v;
     for k=1:numel(v.value)
        w.value{k}=v.value{k}+u;    
     end
     w.s = max(v.s,size(u));   
    else

    if all(size(u)==1)
    value = v.value + u ;
    else
    u=repmat(u(:),1,size(u.value,2));
    value = u.value + v ;
    end
    
    w=v;
    w.s=v.s;
    w.value = value;
    end
    
else
       error('plus not defined') 
end

