function w=multisum(u,k)

if nargin==1
w=reshape(sum(u.value,2),u.s);
else
w = switchmulti(sum(switchmulti(u),k));    
if prod(sizem(u))==1
w = reshape(w.value,w.s);  
end
end