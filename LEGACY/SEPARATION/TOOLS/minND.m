function n = minND(u,varargin)
% function n = minND(u,varargin)
% min o min o ... o min(u)   pour un ND array
s=size(u);
n = u;
for k=1:length(s)
    n = min(n);
end

