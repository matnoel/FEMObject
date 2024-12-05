function s=mrdivide(u,v)

%if isa(u,'MULTIMYDOUBLEND') & isa(v,'MULTIMYDOUBLEND')
%   if u.multidim~=v.multidim | ~all(u.sm==v.sm)
%       error('les dimensions multi doivent correspondre')
%   end
%   s=u;
%   s.value=mrdivide(u.value,v.value);
%   
%elseif ~isa(u,'MULTIMYDOUBLEND')
%   s=v;
%   s.value = mrdivide(u,v.value);
%else
%   s=u;
%   s.value = mrdivide(u.value,v);
%end


if ~israndom(v)
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
    
   error('pas programme') 
end