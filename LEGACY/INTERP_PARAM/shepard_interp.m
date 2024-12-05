function [Jmu,w] = shepard_interp(Jk,muk,mu,p)
% function [Jmu,w] = shepard_interp(Jk,muk,mu)
% Shepard Interpolation
if nargin == 3
    p = 2;
end

n = numel(Jk);
w = zeros(n,1);

for k = 1:n
    w(k) = norm(mu-muk(:,k))^p;
end

m = find(w == 0);

if ~isempty(m)
    Jmu = Jk{m};
else
    w = 1./w;
    Jmu = w(1)*Jk{1};
    for k = 2:n
        Jmu = Jmu+w(k)*Jk{k};
    end
    Jmu = Jmu/sum(w);
end


end