function M = applyfunctiontofaces(M,fun,varargin)
%function M = applyfunctiontofaces(M,fun,varargin)
% fun : function handle
% varargin : arguments de la fonction
% appel de fun(F,varargin{:})
% pour toutes les faces F

for i=1:length(M.facets)
    M.facets{i} = fun(M.facets{i},varargin{:});
end
for i=1:length(M.ridges)
    try
        M.ridges{i} = fun(M.ridges{i},varargin{:});
    end
end
for i=1:length(M.peaks)
    try
        M.peaks{i} = fun(M.peaks{i},varargin{:});
    end
end
