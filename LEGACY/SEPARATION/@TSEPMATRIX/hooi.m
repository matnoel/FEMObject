function [T,err] = hooi(X,r,varargin)
% function T = hooi(X,r,varargin)
% Higher-order orthogonal iteration
% X : TSEPMATRIX tensor
% r : rank of the HOSVD
% T : TSEPMATRIX

gx=gathervectors(X);

[a,err]=hooi(X.alpha,r);

ga=gathervectors(a);
for i=1:X.dim
    ga.F{i}=gx.F{i}*ga.F{i};
end
T=splitvectors(ga);
