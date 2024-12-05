function varargout = plotpeaks(S,varargin)
% function varargout = plotpeaks(S,scanpeaks,varargin)
% scanpeaks : numero des peaks a afficher
% si scanpeaks =0 : toutes les peaks 

Handles = plotfaces(S.peaks,varargin{:});
if nargout>=1
    varargout{1}=Handles;
end