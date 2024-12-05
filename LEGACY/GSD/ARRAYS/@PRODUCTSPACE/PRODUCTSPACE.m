function V = PRODUCTSPACE(dim)

if nargin==0
    dim =[];
    V.dim = dim;
    V = class(V,'PRODUCTSPACE'); 
elseif nargin==1 && isa(dim,'PRODUCTSPACE')
    V = dim;
elseif isa(dim,'double')
    V.dim = dim;
    V = class(V,'PRODUCTSPACE');     
else
    error('constructeur de PRODUCTSPACE attend un double')
end
