function [approx,R,lambda,Xim,indx] = certifiedNonIntrusiveEim(N,fun,snParamFun,Xi,optsEim,K)
% function [approx,R,lambda,Xim,indx] = certifiedNonIntrusiveEim(N,fun,snParamFun,Xi,optsEim,K)

snParamFun0 = snParamFun;
% Extract snapshots of interest such that the submatrix is
% full-rank
m = rank(snParamFun);
[~,~,E] = qr(snParamFun,0);
indxi = E(1:m);
Xim = Xi(:,indxi);
snParamFun = snParamFun(:,indxi); % Snapshots of interest
[Q,R] = qr(snParamFun',0); % Orthogonalization for avoiding
                           % conditioning problems
snParamFun0 = R'\snParamFun0; % Rotation of the original snapshots
snParamFun = Q';

snFun = zeros(N,m);
for k = 1:m
    snFun(:,k) = fun(Xim(:,k));
end
G = snFun/snParamFun;
if nargin == 5
    W = G'*G;
else
    W = G'*K*G;
end
optsEim.metric = @(x) sqrt(x'*W*x);
[~,~,indxi,indx] = eim(snParamFun,optsEim);
lambda = snParamFun(:,indxi);
Xim = Xim(:,indxi);
R = snFun(:,indxi);
lambda = lambda(indx,:) \ snParamFun0(indx,:);
approx = R*lambda;
end
