function s=mrdivide(u,v)


if isa(u,'MULTIMYDOUBLEND') & isa(v,'MULTIMYDOUBLEND')
   s=u;
   if u.multidim~=v.multidim | ~all(u.sm==v.sm)
       error('les dimensions multi doivent correspondre')
   end
   s.value=mrdivide(u.value,v.value);
       
elseif isa(u,'MULTIMYDOUBLEND')
   if isa(v,'MYDOUBLEND') | isa(v,'double')
   if size(v,u.multidim)~=1
      error('la dimension multi doit etre 1 pour le deuxieme argument') 
   end
   s=u;
   s.value = mrdivide(u.value,v);
   else
       error('pas programme')
   end
elseif isa(v,'MULTIMYDOUBLEND')
   if isa(u,'MYDOUBLEND') | isa(u,'double')
   if size(u,v.multidim)~=1
      error('la dimension multi doit etre 1 pour le deuxieme argument') 
   end
   s=v;
   s.value = mrdivide(u,v.value);
   else
       error('pas programme')
   end
              
end

 
