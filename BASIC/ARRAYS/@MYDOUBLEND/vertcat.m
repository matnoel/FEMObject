function u = vertcat(varargin)
% function u = vertcat(varargin)

for i=1:nargin
    varargin{i} = double(varargin{i});
    dim(i) = ndims(varargin{i});
end
iref = find(dim==max(dim));
iref = iref(1);
for i=1:nargin
    varargin{i} = samesizeND(varargin{i},varargin{iref});
end

%%% ATTENTION MATLAB NE GERE PAS ENCORE LES CONCATENATIONS D'EMPTY ARRAY,
%%% on le fait tout seul
if isempty(varargin{1})
    nblig = 0;
    nbcol = size(varargin{1},2);
    for i=1:nargin
        nblig = size(varargin{i},1)+nblig;
        if size(varargin{i},2)~=nbcol
            error('le nombre de colonnes doit correspondre')
        end
    end
    u = zerosND([nblig,nbcol,sizeND(varargin{1})]);
else
    u = MYDOUBLEND(vertcat(varargin{:}));
end