function n=getn(PC,i)

if nargin==1
    n=PC.n;
else
    n=PC.n(i);
end

