function s = getsize(D)
% function s = getsize(D)

x1 = double(getcoord(D.P1));
x2 = double(getcoord(D.P2));

s = x2-x1;
