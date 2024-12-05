function w=mldivide(u,v)

if isa(v,'MULTIMATRIX') && isa(u,'double')
w=v;
    if isa(v.value,'cell')
        for k=1:length(v.value)
          w.value{k} = u\v.value{k};
        end
        
    else
        w.value = u\reshape(v.value,v.s(1),v.s(2)*numelm(v));
        w.value = reshape(w.value,numel(v),numelm(v));
    end
    if all(size(u)==1)
           w.s = size(v);
    elseif all(size(v)==1)
        w.s = size(u);
        else
           w.s = [size(u,1),size(v,2)]; 
    end
        
elseif isa(u,'MULTIMATRIX') && isa(v,'double')

    w=u;
    if isa(u.value,'cell')
        for k=1:length(u.value)
          w.value{k} = u.value{k}\b;
        end  
    else
        for k=1:length(u.value)
          w.value{k} = reshape(u.value(:,k),u.s)\b;
        end     
    end

   if all(size(u)==1)
           w.s = size(v);
   elseif all(size(v)==1)
        w.s = size(u);
   else        
           w.s = [size(u,1),size(v,2)]; 
   end

else    
       error('mldivide not defined') 
end

