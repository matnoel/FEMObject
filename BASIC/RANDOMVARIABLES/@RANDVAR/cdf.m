function P = cdf(u,x)
% function P = cdf(u,x)

P = cdf(u.type,x,u.param{:,2});
P(find(P==1))=1-eps;
P(find(P==0))=eps;
