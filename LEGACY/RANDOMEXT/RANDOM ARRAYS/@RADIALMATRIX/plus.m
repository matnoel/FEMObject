function w=plus(u,v)

if isa(v,'RADIALMATRIX') & isa(u,'RADIALMATRIX')
    w=u;
    w.V = multicat(u.V,v.V);   
    w.L = [u.L ; v.L];
    w.D = sparse([u.D,zeros(u.m,v.m);zeros(v.m,u.m),v.D]);
    w.m = length(w.V);
    
elseif isa(v,'double') 
    w=u;
    w.V = u.V + v;
elseif isa(u,'double')
    w=v;
    w.V = v.V + u;
else
       error('plus not defined')     
end
    
