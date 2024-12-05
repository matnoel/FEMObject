function w=mtimes(u,v)
if isa(v,'logical')
    v=double(v);
end

if isa(u,'logical')
    u=double(u);
end

    
if isa(u,'TIMERADIALMATRIX') & ~isa(v,'TIMERADIALMATRIX')
   if isa(v,'TIMEMATRIX')
       w=u;
     if all(size(v)==1)
         
     w.L = times(w.L,v);
     w.L = TIMEMATRIX(w.L,gettimemodel(v));
          
     else

     w = mtimes(expand(u),v);    
         
     end
     
   elseif isa(v,'double')
     
     w=u;
     w.V = u.V*v;
     try 
         nn = normest(v(:));
     catch
         nn = norm(v);     
     end
     if nn<eps 
     w = TIMERADIALMATRIX(size(w.V),u.TIMEMODEL);    
     end
     
   end
    
elseif isa(v,'TIMERADIALMATRIX') && ~isa(u,'TIMERADIALMATRIX')

    if isa(u,'TIMEMATRIX')
    w=v;
     if all(size(u)==1)
         
     w.L = times(w.L,u);
     w.L = TIMEMATRIX(w.L,gettimemodel(v));
          
     else

     w = mtimes(u,expand(v));    
         
     end


     
    elseif isa(u,'double')
     
     w=v;
     w.V = u*v.V;
     
     if numel(u)<100
     try 
         nn = normest(u(:));
     catch
         nn = norm(u);     
     end
     if nn<eps  
     w = TIMERADIALMATRIX(size(w.V),v.TIMEMODEL);    
     end
     end
     
    end
    
elseif isa(u,'TIMERADIALMATRIX') && isa(v,'TIMERADIALMATRIX')

    
    k=0;
   
    for i=1:u.m
        for j=1:v.m
    k=k+1 ;
    
    V{k} = u.V{i}*v.V{j};
    
    %L{k} = multimtimesold(v.L(j),u.DLmasse{i});
    L{k} = mtimes(u.L(i),v.L(j));
    
        end
    end
    
    w = TIMERADIALMATRIX(V,size(V{1}),L);
   
    
else
       error('mtimes not defined') 
end

