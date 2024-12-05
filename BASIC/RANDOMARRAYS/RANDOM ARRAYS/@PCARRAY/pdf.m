function [P,x] = pdf(u,x)
n=10000;
    us=random(u,n);
if nargin==2
    P=ksdensity(us,x);
else
    [P,x]=ksdensity(us);
end