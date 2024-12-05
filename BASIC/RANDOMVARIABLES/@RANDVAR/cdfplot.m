function varargout = cdfplot(u,varargin)
% function [P,x] = cdfplot(u,varargin)
% trace la cdf de u
% appelle [P,x]=plot(u,'cdf',varargin{:});
varargout=cell(1,nargout);
[varargout{:}]=plot(u,'cdf',varargin{:});