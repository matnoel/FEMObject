function varargout = plotfacets(S,varargin)
% function varargout = plotfacets(S,scanfacets,varargin)
% scanfacets : numero des facets a afficher
% si scanfacets = 0 : toutes les facets 

Handles = plotfaces(S.facets,varargin{:});
if nargout>=1
    varargout{1}=Handles;
end