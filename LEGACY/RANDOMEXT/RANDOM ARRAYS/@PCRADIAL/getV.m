function V=getV(rad,k)

if nargin==1
    V=rad.V;
else
   V=rad.V{k}; 
end