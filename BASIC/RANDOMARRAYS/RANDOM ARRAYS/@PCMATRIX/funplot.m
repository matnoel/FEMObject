function varargout = funplot(upc,varargin)
% function [P,x] = funplot(upc,varargin)
% trace la fonction upc(xi) sur le domaine de definition de la VA xi 
% appelle [P,x]=plot(upc,'fun',varargin{:});
varargout=cell(1,nargout);
[varargout{:}]=plot(upc,'fun',varargin{:});
