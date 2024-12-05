function s=sizem(u,k)
if nargin==1
s=u.sm;
else
s=u.sm(k);
end