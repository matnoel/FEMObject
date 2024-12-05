function w=plus(u,v)

if isa(v,'PCRADIALMATRIX') & isa(u,'PCRADIALMATRIX')
    w=u;
    w.V = multihorzcat(u.V,v.V);   
    w.L = [u.L ; v.L];
    w.D = sparse([u.D,zeros(u.m,v.m);zeros(v.m,u.m),v.D]);
    w.m = length(w.V);
    if ~isempty(u.DLmasse) & ~isempty(v.DLmasse)
    w.DLmasse = multihorzcat(u.DLmasse,v.DLmasse);  
    end
    
elseif isa(u,'PCRADIALMATRIX') & isa(v,'PCMATRIX')
    w=expand(u) + v;
 
elseif isa(v,'PCRADIALMATRIX') & isa(u,'PCMATRIX')
    w=expand(v) + u;

elseif isa(v,'double') 
    w=u+ (v.*ones(1,getPC(u)));
elseif isa(u,'double')
    w=v+(u.*ones(1,getPC(v)));
else
       error('plus not defined')     
end
    
