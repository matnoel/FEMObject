function DN = getDN(elem,xgauss)

DN = get(elem,'DN');
if nargin==2
    DN = DN(xgauss);
end