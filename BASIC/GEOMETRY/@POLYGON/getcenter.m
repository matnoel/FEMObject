function P = getcenter(D)
% function P = getcenter(D)

P = mean(double(squeeze(getcoord(D.P))),1);
P = POINT(P);
