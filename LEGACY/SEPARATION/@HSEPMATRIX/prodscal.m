function w = prodscal(u,v,dim,varargin)

if isa(u,'double')
   w=v;
   for i=1:v.m
       for j=1:length(dim)
      k=dim(j);
      w.F{i,k} = prodscal(u,v.F{i,k});
       end  
   end
elseif isa(v,'double')
   w=u;
   for i=1:u.m
       for j=1:length(dim)
           k=dim(j);
   w.F{i,k} = prodscal(u.F{i,k},v);
       end
   end 
    
else

w = times(u,v);
w = sum(sum(w));

%n = full(reshape([w.F{:}],w.m,w.dim));
w = (sum(prod(n,2).*w.alpha(:)));
end