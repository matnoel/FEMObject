function w=getmatrix(u,k)

if nargin==1 && size(u.value,2)==1
    w=reshape(u.value,u.s);
else
    w=reshape(u.value(:,k),u.s);
end