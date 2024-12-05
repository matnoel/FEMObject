function u = unfreevector(u,varargin)
% function u = unfreevector(u,varargin)

for i=1:length(u.funs)
u.funs{i} = unfreevector(u.funs{i},varargin{:});
end
