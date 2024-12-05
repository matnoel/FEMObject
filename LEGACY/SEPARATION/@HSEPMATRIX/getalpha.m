function m = getalpha(u,i)

if nargin==1
m = u.alpha;
else
m = u.alpha(i);
end
