function w = mtimes(u,v,dim,varargin)

if isa(u,'double')
 if all(size(u)==1)   
   w=v;
   w.alpha = v.alpha*u;
 else  
   w=v;
   for i=1:v.m
   w.F{i,dim} = u.*v.F{i,dim};
   end   
 end
elseif isa(v,'double')
   if all(size(v)==1)   
   w=u;
   w.alpha = u.alpha*v;
   else
   w=u;
   for i=1:u.m
   w.F{i,dim} = u.F{i,dim}.*v;
   end 
   end 
 
    
else
w=u;
w.dim = u.dim;
w.m = u.m*v.m;
w.F = cell(w.m,w.dim);
w.alpha = zeros(1,w.m);

for i=1:u.m
for j=1:v.m
    I=(i-1)*v.m+j;
for k=1:u.dim
w.F{I,k} = u.F{i,k}.*v.F{j,k}; 
end
w.alpha(I) = u.alpha(i)*v.alpha(j); 
end
end


end