function w=mtimes(u,v)

if isa(u,'PCRADIAL') & ~isa(v,'PCRADIAL')
   w=u;
   if isa(v,'PCARRAY') & size(v,1)==1
        for k=1:w.m
         w.L{k} =   v * w.Lmasse{k} ;
         %w.Lmasse{k} = calc_ximasse(getPC(w),w.L{k}); 
        end
    elseif isa(v,'PCARRAY')
       w = double(u.V{1}) * v * u.Lmasse{1} ;
       for k=2:u.m
       w = w + double(u.V{k}) * v * u.Lmasse{k} ;    
       end

   else
        
        for k=1:w.m
         w.V{k} = w.V{k} * v ;   
        end
   
   end
    
elseif isa(v,'PCRADIAL') & ~isa(u,'PCRADIAL')
    w=v;
    if isa(u,'PCARRAY')
        for k=1:w.m
         w.L{k} = u* w.Lmasse{k}  ;
         %w.Lmasse{k} = calc_ximasse(getPC(w),w.L{k}); 
        end
    else
        for k=1:w.m
         w.V{k} = u * w.V{k} ;   
        end
   
    end
    
elseif isa(u,'PCRADIAL') & isa(v,'PCRADIAL')
    k=0;
    
    for i=1:u.m
        for j=1:v.m
    k=k+1 ;
    V{k} = u.V{i}*v.V{j};
    L{k} = v.L{j}*u.Lmasse{i};
        end
    end
    w = PCRADIAL(V,L)
else
       error('mtimes not defined') 
end

