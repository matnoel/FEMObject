function connecvalue = getconnecvalue(elem,value,node)

if nargin==3
connec = calc_conneclocal(elem,node);
else
connec = getconnec(elem);
end
connecvalue = double(value);
if size(connecvalue,2)==1
connecvalue = connecvalue(:);    
end
connecvalue = reshape(connecvalue(connec,:),[size(connec),size(connecvalue,2)]);


