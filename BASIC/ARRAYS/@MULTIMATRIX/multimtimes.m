function w=multimtimes(u,v)
%warning('nouvelle fonction')

if (isa(u,'MULTIMATRIX') && isa(u.value,'cell') ) || (isa(v,'MULTIMATRIX') && isa(v.value,'cell'))

    
    if isa(v,'double')
    
    for i=1:size(u.value,1)
    for j=1:size(v,2)
    value{i,j} = u.value{i,1}*v(1,j);
    for k=2:size(v,1)
    value{i,j} = value{i,j} + u.value{i,k}*v(k,j);    
    end
    end
    end
    
    w = MULTIMATRIX(value,size(u),[size(u.value,1),size(v,2)]);
    elseif isa(u,'double')
        
       
    
    for i=1:size(u,1)
    for j=1:size(v.value,2)
    value{i,j} = sparsemtimes(u(i,1),v.value{1,j});
    for k=2:size(v.value,1)
    value{i,j} = value{i,j} + sparsemtimes(u(i,k),v.value{k,j});    
    end
    end
    end
    w = MULTIMATRIX(value,size(v),[size(u,1),size(v.value,2)]);
        
    elseif isa(v,'MULTIMATRIX') && isa(v.value,'double') 
        
        w = multimtimes(u,mat2cell(v));
    elseif isa(u,'MULTIMATRIX') && isa(u.value,'double')
        
        w = multimtimes(mat2cecll(u),v);        
    else
        if ~all(u.sm(2)==v.sm(1))
            error('pas les meme multidim')
        end
    w = v ;
    
    for i=1:size(u.value,1)
    for j=1:size(v.value,2)
    value{i,j} = sparsemtimes(u.value{i,1},v.value{1,j});
    for k=2:size(v.value,1)
    value{i,j} = value{i,j} + sparsemtimes(u.value{i,k},v.value{k,j});    
    end
    end
    end
    w = MULTIMATRIX(value,size(v),[size(u.value,1),size(v.value,2)]);

    end
    
    
else  
    
if isa(u,'MULTIMATRIX')
    u = switchmulti(u);
end
if isa(v,'MULTIMATRIX')
    v = switchmulti(v);
end

w=switchmulti(mtimes(u,v));
end
%if isa(u,'MULTIMATRIX') & isa(v,'double')
%    w=u;
%    w.value = u.value * v ;
%    w.sm = [size(w.value,2),1];
%end
