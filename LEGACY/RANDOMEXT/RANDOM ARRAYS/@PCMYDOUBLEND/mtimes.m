function s=mtimes(u,v)

if ~israndom(u)
   if isa(u,'double') | isa(u,'MYDOUBLEND')
      s=v;
      if size(u,v.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = mtimes(u,v.V);
    else
        error('pas programme')
    end
elseif ~israndom(v)
    
    if isa(v,'double') | isa(v,'MYDOUBLEND')
      s=u;
      if size(v,u.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = mtimes(u.V,v);
    else
        error('pas programme')
    end
else
error('pas programme')
    if isa(u,'PCMYDOUBLEND')
    if isradial(u)
        up=PCRADIALMATRIX(u);
    else
        up=PCMATRIX(u);
    end
    up=reshape(up,prod(size(up)),1);
    else
    up =u;
    end
    
    if isa(v,'PCMYDOUBLEND')
    if isradial(v)
        vp=PCRADIALMATRIX(v);
    else
        vp=PCMATRIX(v);
    end
    vp=reshape(vp,prod(size(vp)),1);
    else
    vp =v;
    end
    
    s = mtimes(up,vp);
    
    
end
     
     