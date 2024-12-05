function x = mtimes(a,b)

if isa(a,'PCTPMATRIX') && isa(b,'double')
    x = a;
    x.phi0 = mtimes(a.phi0,b);
elseif isa(b,'PCTPMATRIX') && isa(a,'double')
    x = b;
    x.phi0 = mtimes(a,b.phi0);
elseif isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIX')
    
    if ~isximasse(a) && ~isximasse(b)
    warning('pourrait etre optimise en precalculant ximasse')  
    a = calc_ximasse(a);
    end
    x=a;
    
    for i=1:getnbdim(a)
    isranda = isranddim(a,i);
    israndb = isranddim(b,i);
    if isranda && israndb
    if isximasse(a)    
    x.phi{i} = get_ximasse(a,i)*b.phi{i};  
    else
    x.phi{i} = get_ximasse(b,i)*a.phi{i};    
    end
    else       
    x.phi{i} = a.phi{i}*b.phi{i};      
    end                   
    end
    
    x.phi0 = a.phi0*b.phi0;
    x.ximasse={};
    x.isranddim = isranddim(a) | isranddim(b);
else
    error('pas prevu')
    
end
