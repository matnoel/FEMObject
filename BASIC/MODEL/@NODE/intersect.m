function [N,numnode1,numnode2] = intersect(N1,N2)
% function [N,numnode1,numnode2] = intersect(N1,N2)

[P,rep1,rep2]=intersect(POINT(N1),POINT(N2));
numnode1=getnumber(N1,rep1);
numnode2=getnumber(N2,rep2);
N=getnode(N1,numnode1);