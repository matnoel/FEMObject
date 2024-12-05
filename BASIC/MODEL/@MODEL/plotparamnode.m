function varargout = plotparamnode(M,varargin)
% function [Handles,leg] = plotparamnode(M,choix,varargin)
% choix : parametres à afficher
% Handles : handles vers les patchs d'elements
% leg : cellules contenant les legendes de chaque handle

plot(M);
varargout = cell(1,nargout);
[varargout{:}] = plotparamnode(getnode(M),varargin{:});
dilateaxis(1.1)

