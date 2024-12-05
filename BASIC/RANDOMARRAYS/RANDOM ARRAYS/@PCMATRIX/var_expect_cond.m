function s = var_expect_cond(u,dim)
% function s = var_expect_cond(u,dim)
% compute the variance of the conditional expectation of u with respect to
% variables indicated in vector dim

PC = getPC(u);
moy = expect(u);

if nargin==1
    dim = 1:getM(PC);
end

s = zeros(length(moy),length(dim));
if getM(PC)==1
    var = expecttimes(u,u)-moy.^2;
    s(:,1) = var;
else
    for k=1:length(dim)
        i = dim(k);
        mi = expectnodim(i,u);
        vi = expecttimes(mi,mi)-moy.^2;
        s(:,k) = vi;
    end
end
