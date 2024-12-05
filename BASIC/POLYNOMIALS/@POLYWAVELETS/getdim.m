function n = getdim(h,p)
if nargin==1
    p=getparam(h,'p');
end
n=(p+1)*(2^(getparam(h,'n')+1));


