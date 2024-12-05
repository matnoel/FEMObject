function an = norm(a,A)
% function an = norm(a,A)

b = a;
pcdim = a.pcdim;
a = double(a);
na = ndims(a);
a = permute(a,[pcdim,setdiff(1:na,pcdim)]);
sa = size(a);
if nargin==2
    b = a*A;
else
    b = a;
end

an = sqrt(sum(sum(a(:,:).*b(:,:),2),1));
