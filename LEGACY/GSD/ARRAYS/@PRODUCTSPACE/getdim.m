function dim = getdim(V,j)

if nargin==1
dim = V.dim;
else
dim = V.dim(j);    
end