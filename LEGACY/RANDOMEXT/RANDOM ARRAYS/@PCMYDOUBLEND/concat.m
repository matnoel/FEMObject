function u = concat(u,v,k)

if isa(u,'PCMYDOUBLEND') && isa(v,'PCMYDOUBLEND')

    if isradial(u) && ~isradial(v)
      u=expand(u);
    end
    if isradial(v) && ~isradial(u)
      v=expand(v);
    end
    u.V = concat(u.V,v.V,k);
     
elseif isa(u,'PCMYDOUBLEND') && isa(v,'MYDOUBLEND')
    
  if isradial(u)  
   u.V = concat(u.V,v,k);   
  else
   v = v*one(getPC(u));   
   u = concat(u,v,k)
  end
    
elseif isa(v,'PCMYDOUBLEND') && isa(u,'MYDOUBLEND')
    
  if isradial(v)      
   v.V = concat(u,v.V,k);
   u=v;
  else
   u = u*one(getPC(u));   
   u = concat(u,v,k);
  end
    
else
    
    error('pas prevu')
    
end