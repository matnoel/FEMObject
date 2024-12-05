function w = sparsetimes(u,v)

if issparse(u) & all(size(v)==1)
    w = spfun(@(x) x.*v,u); 
elseif issparse(v) & all(size(u)==1)
    w = spfun(@(x) u.*x,v);     
else
    w = times(u,v); 
end
