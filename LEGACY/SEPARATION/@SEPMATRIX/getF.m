function F = getF(u,m,k)
% function F = getF(u,m,k)
% m: index of rank-one functions
% k: dimension
if nargin==1
F = u.F;
elseif nargin==2 
F = u.F(m,:);
elseif nargin==3 && isempty(m)
F = u.F(:,k);
elseif nargin==3
F = u.F{m,k};
end

