function [u,varargout] = dsolve(M,A,b,N,varargin)

varargout=cell(1,nargout-1);
[u,varargout{:}] = dsolve_alterne(M,A,b,N,varargin{:});