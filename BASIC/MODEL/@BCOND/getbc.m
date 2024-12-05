function b = getbc(B,i)

if nargin==2
    b = B.BC{i};
else
    b = B.BC;
end