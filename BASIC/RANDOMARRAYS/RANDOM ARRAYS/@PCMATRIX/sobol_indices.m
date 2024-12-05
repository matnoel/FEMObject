function s = sobol_indices(u,dim)
% function s = sobol_indices(u,dim)
% compute the first order (partial) sensitivity indices of u with respect
% to variables indicated in vector dim

PC = getPC(u);
moy = expect(u);

if nargin==1
    dim = 1:getM(PC);
end

s = zeros(length(moy),length(dim));
if getM(PC)==1
    s(:,1) = 1;
else
    var = expecttimes(u,u)-moy.^2;
    for k=1:length(dim)
        i = dim(k);
        mi = expectnodim(i,u);
        vi = expecttimes(mi,mi)-moy.^2;
        s(:,k) = vi./var;
    end
end
