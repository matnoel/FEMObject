function S = initializeBN(S,varargin)

for p=1:S.nbgroupelem
S.groupelem{p} = initializeBN(S.groupelem{p},S.node,varargin{:});
end