function metric=getmetric(SM,dim)
if nargin==1
    metric=SM.metric;
elseif nargin==2
    metric=SM.F{dim}.metric;
end