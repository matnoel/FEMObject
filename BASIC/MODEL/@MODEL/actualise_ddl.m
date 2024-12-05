function M = actualise_ddl(M,varargin)
% function M = actualise_ddl(M,varargin)

for p=1:length(M.groupelem)
    M.groupelem{p} = actualise_ddl(M.groupelem{p},varargin{:});
end

