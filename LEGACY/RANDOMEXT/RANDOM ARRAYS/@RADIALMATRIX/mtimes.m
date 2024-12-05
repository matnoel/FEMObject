function w=mtimes(u,v)

if isa(u,'PCRADIAL') & ~isa(v,'PCRADIAL')
   if isa(v,'double') 
   w=u;
   w.V = u.V * v ;
   elseif isa(v,class(u.L))
   w=u;
   w.L=u.L * v ;
   end
    
elseif isa(v,'PCRADIAL') & ~isa(u,'PCRADIAL')
   if isa(u,'double') 
   w=v;
   w.V = u * v.V ;
   elseif isa(u,class(vu.L))
   w=v;
   w.L= u * v.L;
    end
   
elseif isa(u,'PCRADIAL') & isa(v,'PCRADIAL')
    k=0;
    for i=1:u.m
    for j=1:v.m
    k=k+1 ;
    V{k} = u.V{i}*v.V{j};
    L{k} = u.L{i}*v.L{j};
        end
    end
    w = RADIALMATRIX(V,L,size(V{1}))
else
       error('mtimes not defined') 
end

