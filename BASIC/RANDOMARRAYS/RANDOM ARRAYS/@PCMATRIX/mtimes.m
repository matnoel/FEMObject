function w=mtimes(u,v)
if isa(u,'logical')
    u=double(u);
elseif isa(v,'logical')
    v=double(v);
end
if isa(u,'PCMATRIX') && isa(v,'PCMATRIX')

    PCu = getPC(u);
    PCv = getPC(v);
    if polycmp(PCv,PCu)    
    masse= MULTIMATRIX(getmasse(PCu));
    if isempty(masse)
    masse= MULTIMATRIX(getmasse(PCv));    
    value = mtimes(u.MULTIMATRIX,v.MULTIMATRIX,masse,1);
    else
    value = mtimes(u.MULTIMATRIX,v.MULTIMATRIX,masse,2);                
    end
    elseif isin(PCv,PCu)
    v = project(v,PCu);
    %PCu = calc_masse(PCu,PCv);
    PCu = calc_masse(PCu);
    masse= MULTIMATRIX(getmasse(PCu));
    value = mtimes(u.MULTIMATRIX,v.MULTIMATRIX,masse,1);
    u.MULTIMATRIX = value;
    u = project(u,PCv);
    value = u.MULTIMATRIX;
    else
    PCv = calc_masse(PCv,PCu);    
    masse= MULTIMATRIX(getmasse(PCv));
    value = mtimes(u.MULTIMATRIX,v.MULTIMATRIX,masse,2);        
    end
    w=u;
    w.MULTIMATRIX = value;
    
    
elseif isa(u,'PCMATRIX') && isa(v,'PCARRAY')
   w = u .* PCMATRIX(v) ;
elseif isa(u,'PCARRAY') && isa(v,'PCMATRIX')
   w = PCMATRIX(u) * v ;
elseif isa(u,'PCMATRIX') && isa(v,'double')
    if all(size(u)==1) && ~all(size(v)==1)
    if all(size(v)>1)
    w=PCRADIALMATRIX({v},size(v),u);
    else
    w=PCRADIALMATRIX(v(:),size(v),u);    
    end
    else
    w=u;
    w.MULTIMATRIX=u.MULTIMATRIX * v ;
    end
elseif isa(v,'PCMATRIX') && isa(u,'double')
    if all(size(v)==1) && ~all(size(u)==1)
    if all(size(u)>1) 
    w=PCRADIALMATRIX({u},size(u),v);    
    else
    w=PCRADIALMATRIX(u(:),size(u),v);
    end
    else
    w=v;
    w.MULTIMATRIX=u*v.MULTIMATRIX  ;
    end
elseif isa(u,'PCMATRIX') && isa(v,'MYDOUBLEND')
    if all(size(u)==1)
    w=PCMYDOUBLEND(u,5)*v;
    else
    error('not defined')
    end
elseif isa(v,'PCMATRIX') && isa(u,'MYDOUBLEND')
    if all(size(v)==1)
    w=u*PCMYDOUBLEND(v,5);
    else
    error('not defined')
    end
elseif isa(u,'cell') && isa(v,'PCMATRIX') 
    if all(size(v)==1)
    w=PCRADIALMATRIX(u,size(u{1}),v);
    else
    error('not defined')
    end

elseif isa(u,'PCMATRIX') && isa(v,'cell')
    if all(size(u)==1)
    w=PCRADIALMATRIX(v,size(v{1}),u);
    else
    error('not defined')
    end
    
else
       error('mtimes not defined') 
end

