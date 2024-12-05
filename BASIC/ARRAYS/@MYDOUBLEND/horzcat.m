function u = horzcat(varargin)
% function u = horzcat(varargin)

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
    nbcol = 0;
    nblig = size(varargin{1},1);
    for i=1:nargin
        nbcol = size(varargin{i},2)+nbcol;
        if size(varargin{i},1)~=nblig
            error('le nombre de lignes doit correspondre')
        end
    end
    u = zerosND([nblig,nbcol,sizeND(varargin{1})]);
else
    u = MYDOUBLEND(horzcat(varargin{:}));
end
