function s=mrdivide(u,v)

if ~israndom(u)
   if isa(u,'double') | isa(u,'MYDOUBLEND')
      s=v;
      if size(u,v.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = mrdivide(u,v.V);
    else
        error('pas programme')
    end
elseif ~israndom(v)
    
    if isa(v,'double') | isa(v,'MYDOUBLEND')
      s=u;
      if size(v,u.stodim)>1
          error('les dimensions ne correspondent pas')
      end
      s.V = mrdivide(u.V,v);
    else
        error('pas programme')
    end
else
error('multiplication non possible')
end
     
     