function a = quantile(apc,q,n)

a = random(apc,n);

a = quantile(double(a)',q)';
if length(q)>1
a = MULTIMATRIX(a,size(apc),[length(q),1]);
else
a = reshape(a,size(apc));    
end