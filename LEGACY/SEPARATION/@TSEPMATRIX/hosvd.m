function [T,err] = hosvd(X,r)
% function [T,err] = hosvd(X,r)
% X : TSEPMATRIX
% r : rank of the HOSVD
% T : TSEPMATRIX

gx=gathervectors(X);

[a,err]=hosvd(X.alpha,r);
ga=gathervectors(a);
for i=1:X.dim
    ga.F{i}=gx.F{i}*ga.F{i};
end
T=splitvectors(ga);
