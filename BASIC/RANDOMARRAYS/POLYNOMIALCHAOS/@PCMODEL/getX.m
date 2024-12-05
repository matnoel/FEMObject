function X = getX(PC,i)

if nargin==1
X=PC.X;
else
X = PC.X(i);    
end
