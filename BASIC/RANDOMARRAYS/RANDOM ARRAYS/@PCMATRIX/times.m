function w=times(u,v)

if isa(u,'PCMATRIX') && isa(v,'PCMATRIX')
    PCu = getPC(u);
    PCv = getPC(v);
    if polycmp(PCv,PCu)    
    masse= MULTIMATRIX(getmasse(PCu));
    if isempty(masse)
    masse= MULTIMATRIX(getmasse(PCv));    
    value = times(u.MULTIMATRIX,v.MULTIMATRIX,masse,1);
    else
    value = times(u.MULTIMATRIX,v.MULTIMATRIX,masse,2);                
    end
    elseif isin(PCv,PCu)
    masse= MULTIMATRIX(getmasse(PCu));
    value = times(u.MULTIMATRIX,v.MULTIMATRIX,masse,1);    
    else
    masse= MULTIMATRIX(getmasse(PCu));
    value = times(u.MULTIMATRIX,v.MULTIMATRIX,masse,2);        
    end
    w=u;
    w.MULTIMATRIX = value;
    
    
elseif isa(u,'PCMATRIX') && isa(v,'PCARRAY')
   w = u .* PCMATRIX(v) ;
elseif isa(u,'PCARRAY') && isa(v,'PCMATRIX')
   w = PCMATRIX(u) .* v ;
elseif isa(u,'PCMATRIX') && isa(v,'double')
    if all(size(u)==1) && ~all(size(v)==1)
    if all(size(v)>1)
    w=PCRADIALMATRIX({v},size(v),u);
    else
    w=PCRADIALMATRIX(v(:),size(v),u);    
    end
    else
    w=u;
    w.MULTIMATRIX=u.MULTIMATRIX .* v ;
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
    w.MULTIMATRIX=u .* v.MULTIMATRIX  ;
    end

else
    
       error('times not defined') 
end

