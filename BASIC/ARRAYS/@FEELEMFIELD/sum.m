function u = sum(u,varargin)
% function u = sum(u,varargin)

for k=1:length(u.value)
    u.value{k} = sum(u.value{k},varargin{:});
end