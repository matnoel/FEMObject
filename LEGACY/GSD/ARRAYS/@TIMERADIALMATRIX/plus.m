function w=plus(u,v)

if isa(v,'TIMERADIALMATRIX') && isa(u,'TIMERADIALMATRIX')
    w=u;

    w.V = multihorzcat(u.V,v.V);   
    w.L = [u.L ; v.L];
    w.D = sparse([u.D,zeros(u.m,v.m);zeros(v.m,u.m),v.D]);
    w.m = length(w.V);
    
    
elseif isa(u,'TIMERADIALMATRIX') && isa(v,'TIMEMATRIX')
    w=expand(u) + v;
 
elseif isa(v,'TIMERADIALMATRIX') && isa(u,'TIMEMATRIX')
    w=expand(v) + u;

elseif isa(v,'double') 
    if normest(v)==0
    w=u;
    else
    w=u+ TIMERADIALMATRIX(v,size(v),ones(1,gettimemodel(u)));
    end
elseif isa(u,'double')
    if normest(u)==0
    w=v;
    else
    w=v+TIMERADIALMATRIX(u,size(u),ones(1,gettimemodel(v)));
    end
    
else
       error('plus not defined')     
end
    
