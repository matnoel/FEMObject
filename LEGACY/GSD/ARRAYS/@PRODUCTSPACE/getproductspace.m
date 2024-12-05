function V = getproductspace(V,dim)

if ~isa(V,'PRODUCTSPACE')
    V = PRODUCTSPACE(getdim(V));
elseif nargin>=2
    
if ~all(ismember(dim,V.dim))
    error('les dimensions indiquees ne sont pas dans le PRODUCTSPACE')
end

V.dim = dim;
end


