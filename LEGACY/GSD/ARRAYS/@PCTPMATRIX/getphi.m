function phi = getphi(u,i)

if nargin==1
phi=u.phi;
elseif length(i)==1
phi=u.phi{i};    
else
phi=u.phi(i);      
end
