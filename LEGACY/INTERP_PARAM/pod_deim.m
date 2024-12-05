function [u,ind] = pod_deim(u,tol)
% function [u,ind] = pod_deim(u,tol)

% POD
[u,s,~] = svd(u,0);
s = diag(s);
e = 1-sqrt(cumsum(s.^2)/sum(s.^2));
m = find(e < tol,1);
u = u(:,1:m);

% DEIM
n = size(u,1);
ind = zeros(m,1);
[~,ind(1)] = max(abs(u(:,1)));

for l = 2:m
    PtU = u(:,1:l-1);
    PtU = PtU(ind(1:l-1),:);
    c = PtU\u(ind(1:l-1),l);
    r = u(:,l)-u(:,1:l-1)*c;
    
    [~,ind(l)] = max(abs(r));
end

end