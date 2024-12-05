function w = mldivide(u,v,dim,varargin)

if nargin==3 && isa(u,'double') 
w=v;
for j=1:v.m
w.F{j,dim} = u\v.F{j,dim}; 
end

elseif isa(u,'double')
 if all(size(u)==1)   
   w=v;
   w.alpha = u\v.alpha;
 elseif nargin==2 
   error('pas prevu, inverser l''operation')    
 else
   w=v;
   for i=1:v.m
      w.F{i,dim} = u\v.F{i,dim};
   end
 end
elseif isa(v,'double')  
    error('pas prevu, inverser l''operation')
else
w=u;
w.dim = u.dim;
if u.m>1
    error('first argument must be rank one')
end
w.m = v.m;
w.F = cell(w.m,w.dim);
w.alpha = zeros(1,w.m);
if nargin==2
   dim = 1:u.dim; 
end
for j=1:v.m
for kkk=1:length(dim)
k=dim(kkk);
w.F{j,k} = u.F{1,k}\v.F{j,k}; 
end
w.alpha(j) = u.alpha\v.alpha(j); 
end


end