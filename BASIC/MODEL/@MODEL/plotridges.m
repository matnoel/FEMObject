function varargout = plotridges(S,varargin)
% function varargout = plotridges(S,scanridges,varargin)
% scanridges : numero des ridges a afficher
% si scanridges = 0 : toutes les ridges 

Handles = plotfaces(S.ridges,varargin{:});
if nargout>=1
    varargout{1}=Handles;
end