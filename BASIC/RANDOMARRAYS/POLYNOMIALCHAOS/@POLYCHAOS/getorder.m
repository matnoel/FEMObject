function p=getorder(PC,i)

if nargin==1
    p=PC.p;
else
    p=PC.p(i);
end
