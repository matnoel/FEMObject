function s=plus(u,v)


if isa(u,'PCMYDOUBLEND') && isa(v,'PCMYDOUBLEND')
    if u.stodim~=v.stodim
        error('les dimensions stochastiques ne correspondent pas')
    end
    
    if ~isradial(u) && ~isradial(v)
        s = u;
        s.V = u.V + v.V ; 
    elseif isradial(u) && isradial(v)
        s = u;
        s.L = vertcat(u.L,v.L);
        s.V = concat(u.V,v.V,u.stodim);
    elseif isradial(u)
        if numel(v)==1
        s = plus(u,convertradial(v));
        else
        s = plus(expand(u),v);    
        end
    elseif isradial(v)
        if numel(u)==1
        s = plus(convertradial(u),v);
        else
        s = plus(u,expand(v));    
        end
    end
    

elseif isa(u,'PCMYDOUBLEND')
    if isa(v,'double')
      v = times(v,one(getPC(u)));
    end
    if isa(v,'PCMATRIX') || isa(v,'PCRADIALMATRIX')
     s = plus(u,PCMYDOUBLEND(v,u.stodim));
    else
        error('pas programme')
    end
elseif isa(v,'PCMYDOUBLEND')
    if isa(u,'double')
      u = times(u,one(getPC(v)));
    end
    if isa(u,'PCMATRIX') || isa(u,'PCRADIALMATRIX')
     s = plus(PCMYDOUBLEND(u,v.stodim),v);
        
    else
        error('pas programme')
    end
else
    error('addition non possible')
end
     
     