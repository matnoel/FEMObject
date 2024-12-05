function s=size(u,k)
if nargin==1
s=u.s;
else
s=u.s(k);
end