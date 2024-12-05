function c = getfun(c,i)

if nargin==1
    c = c.tensorfuns;
else
    c = c.tensorfuns(i);
    if length(i)==1
        c=c{1};
    end
end
    