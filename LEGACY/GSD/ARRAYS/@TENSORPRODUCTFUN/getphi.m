function phi = getphi(a,i)

if nargin==1
    phi=a.phi;
elseif length(i)>1
    phi = a.phi(i);
else
    phi = a.phi{i};
end


