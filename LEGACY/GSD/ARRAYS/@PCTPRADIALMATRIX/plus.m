function w=plus(u,v)

if isa(v,'PCTPRADIALMATRIX') && isa(u,'PCTPRADIALMATRIX')
    w=u;
    
    w.V = multicat(u.V,v.V);   
    w.L = [u.L ; v.L];
    w.D = sparse([u.D,zeros(u.m,v.m);zeros(v.m,u.m),v.D]);
    w.m = length(w.V);
    if ~isempty(u.DLmasse) && ~isempty(v.DLmasse)
    w.DLmasse = multicat(u.DLmasse,v.DLmasse);  
    end
    
    
elseif isa(u,'PCTPRADIALMATRIX') && (isa(v,'PCTPMATRIX') || isa(v,'PCTPMATRIXSUM')) 
    if iscell(v) || iscell(u)
        v=mat2cell(v);
        u=mat2cell(u);
    end
    w=expand(u) + v;
 
elseif isa(v,'PCTPRADIALMATRIX') && (isa(u,'PCTPMATRIX') || isa(u,'PCTPMATRIXSUM')) 
    if iscell(v) || iscell(u)
        v=mat2cell(v);
        u=mat2cell(u);
    end
    w=expand(v) + u;

elseif isa(v,'double') 
    if normest(v)==0
    w=u;
    else
    w=u+ PCTPRADIALMATRIX(v,size(v),ones(1,getPC(u)));
    end
elseif isa(u,'double')
    if normest(u)==0
    w=v;
    else
    w=v+PCTPRADIALMATRIX(u,size(u),ones(1,getPC(v)));
    end
    
else
       error('plus not defined')     
end
    
