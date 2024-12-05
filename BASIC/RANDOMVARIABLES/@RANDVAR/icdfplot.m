function varargout = icdfplot(u,varargin)
% function [P,x] = icdfplot(u,varargin)
% trace la icdf de u
% appelle [P,x]=plot(u,'icdf',varargin{:});
varargout=cell(1,nargout);
[varargout{:}]=plot(u,'icdf',varargin{:});