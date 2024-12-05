function w=mtimes(u,v)
if isa(v,'logical')
    v=double(v);
end

if isa(u,'logical')
    u=double(u);
end


    
if isa(u,'PCTPRADIALMATRIX') && ~isa(v,'PCTPRADIALMATRIX')
   if isa(v,'PCTPMATRIX') || isa(v,'PCTPMATRIXSUM')
     error('pas programme')
     
   elseif isa(v,'double')     
     w=u;
     w.V = u.V*v;
   end
    
elseif isa(v,'PCTPRADIALMATRIX') && ~isa(u,'PCTPRADIALMATRIX')

    if isa(u,'PCTPMATRIX') || isa(u,'PCTPMATRIXSUM')
    error('pas programme')
    
     
    elseif isa(u,'double')
     
     w=v;
     w.V = u*v.V;
     
    end
    
elseif isa(u,'PCTPRADIALMATRIX') && isa(v,'PCTPRADIALMATRIX')
 k=0;
   
    V = cell(1,u.m*v.m);
    L = cell(1,u.m*v.m);
    for i=1:u.m
        for j=1:v.m
    k=k+1 ;    
    V{k} = u.V{i}*v.V{j}; 
    L{k} = mtimes(u.L{i},v.L{j});
        end
    end
    
    w = PCTPRADIALMATRIX(V,size(V{1}),L);
    if isdouble(u) && isdouble(v)
    w = cell2mat(w);   
    end
    
else
       error('mtimes not defined') 
end

