function s=times(u,v)

if ~israndom(u)
   if isa(u,'double') | isa(u,'MYDOUBLEND')
      s=v;
      if size(u,v.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = times(u,v.V);
    else
        error('pas programme')
    end
elseif ~israndom(v)
    
    if isa(v,'double') | isa(v,'MYDOUBLEND')
      s=u;
      if size(v,u.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = times(u.V,v);
    else
        error('pas programme')
    end
else
    if isa(u,'PCMYDOUBLEND') &  isa(v,'PCMYDOUBLEND') 
    if ~all(size(u)==1) && ~all(size(v)==1)
        if length(size(u))~=length(size(v)) || ~all(size(u)==size(v))
        [u,v]=samesize(u,v); 
        end
    end
    end
   
    
    if isa(u,'PCMYDOUBLEND')
    ku = u.stodim;
    su = size(u);
    if isradial(u)
        up=PCRADIALMATRIX(u);
    else
        up=PCMATRIX(u);
    end
    up=reshape(up,prod(size(up)),1);
    else
    up = u ;
    ku = v.stodim;
    end
    
    if isa(v,'PCMYDOUBLEND')
    kv = v.stodim;
    sv = size(v);
    if isradial(v)
        vp=PCRADIALMATRIX(v);
    else
        vp=PCMATRIX(v);
    end
    vp=reshape(vp,prod(size(vp)),1);
    else
    vp = v ;
    kv = u.stodim;
    end
    
    if ku~=kv
            error('les dimensions stochastiques doivent correspondre')
    else
        k=ku;
    end
  
    if prod(su)==1
    ss = sv;
    s = mtimes(up,vp);
    elseif prod(sv)==1
    ss = su;    
    s = mtimes(up,vp);
    else
        if length(su)~=length(sv) || any(su~=sv)
        error('le nombre d''elements doit etre le meme')
        else
        ss=su;
        end
   
    s = times(up,vp);
    end
   
    s = reshape(s,ss);
   
    s = PCMYDOUBLEND(s,ku);
    
     
end
     
     