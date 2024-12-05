function d = pdfdistance(a,b,varargin)
% function d = pdfdistance(a,b,varargin)

m = getcharin('npts',varargin,100);
n = getcharin('nbs',varargin,1e5);

if ~israndom(a)
    as = a ; 
    bs = random(b,n,1);   
elseif ~israndom(b)
    bs = b ; 
    as = random(a,n,1);
end
ax = getcharin('ax',varargin,[min(min(as),min(bs)),max(max(as),max(bs))]);       

x=linspace(ax(1),ax(2),m);

if ischarin('ks',varargin)
    [pa,xm] = ksdensity(as,x);
    [pb,xm] = ksdensity(bs,x);
else
    [pa,xm] = pdfsample(as,x);
    [pb,xm] = pdfsample(bs,x);   
end
pa = max(pa,eps);
pb = max(pb,eps);

d = (kldiv(xm,pa,(pb+pa)/2)+kldiv(xm,pb,(pa+pb)/2))/2;


