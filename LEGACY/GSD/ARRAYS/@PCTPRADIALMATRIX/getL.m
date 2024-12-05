function L=getL(rad,k)
if nargin==1
L=rad.L;
else
L=rad.L{k};    
end