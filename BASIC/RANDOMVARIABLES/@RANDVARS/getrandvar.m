function rv = getrandvar(rvs,k)

if nargin==1
rv = rvs.RV;
else
rv = rvs.RV{k};    
end