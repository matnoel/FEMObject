function metric=calc_metric(h,p,varargin)
%function metric=calc_metric(h,p,varargin)
% calcul de la metric matrix
% metric : array (p+1)-by-(p+1) 
% E(a * b) = a'*metric*b  ou a et b sont les vecteurs
% representant les coeff de a et b sur la base (h0, ... , hp) 

if isorthonormal(h)
metric = speye(getdim(h,p),getdim(h,p));
else
n = getdim(h,p);
metric=sparse(n,n);
for i=0:n-1
for j=i:n-1
metric(i+1,j+1)=moment(h,[i,j]);
metric(j+1,i+1)=metric(i+1,j+1);
end
end
end