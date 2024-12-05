function n = maxND(u,varargin)
% function n = maxND(u,varargin)
% max o max o ... o max(u)   pour un ND array
s=size(u);
n = u;
for k=1:length(s)
    n = max(n);
end

