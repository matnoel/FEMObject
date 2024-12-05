function s = sensitivity_indices_max_expect(u,dim)
% function s = sensitivity_indices_max_expect(u,dim)
% compute the first order (partial) sensitivity indices of u with respect
% to variables indicated in vector dim

PC = getPC(u);
moy = expect(u);

if nargin==1
    dim = 1:getnbgroups(PC);
end

s = zeros(length(moy),length(dim));
if getM(PC)==1
    var = expecttimes(u,u)-moy.^2;
    s(:,1) = var./max(var);
else
    for k=1:length(dim)
        i = dim(k);
        mi = expectnodim(i,u);
        vi = expecttimes(mi,mi)-moy.^2;
        s(:,k) = vi./max(vi);
    end
end
