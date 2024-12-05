function w=plus(u,v)

if isa(v,'PCRADIAL') & isa(u,'PCRADIAL')
    w=u;
    w.V = [u.V v.V];   
    w.L = [u.L v.L];
    w.D = sparse([u.D,zeros(u.m,v.m);zeros(v.m,u.m),v.D]);
    w.m = length(w.V);
    
elseif isa(u,'PCRADIAL') & isa(v,'PCARRAY')
    w=expand(u) + v;
 
elseif isa(v,'PCRADIAL') & isa(u,'PCARRAY')
    w=expand(v) + u;

elseif isa(v,'double') 
    w=u+ (v.*one(getPC(u)));
elseif isa(u,'double')
    w=v+(u.*one(getPC(v)));
else
       error('plus not defined')     
end
    
