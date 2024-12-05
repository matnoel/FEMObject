function e = isenrichtip(c,k)
% function e = isenrichtip(c,k)

if nargin==1
    k=1:getnbtip(c);
end
e = zeros(1,length(k));
for i=1:length(k)    
e(i) = isenrich(c.LEVELSETS{1+k(i)});
end



