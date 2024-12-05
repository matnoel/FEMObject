function varargout = cdfsampleplot(us,varargin)
% function [P,x] = cdfsampleplot(us,varargin)
% trace la cdf d'un ensemble de realisations
% appelle [P,x]=sampleplot(us,'cdf',varargin{:});
%
% See also sampleplot, pdfsampleplot

varargout=cell(1,nargout);
[varargout{:}]=sampleplot(us,'cdf',varargin{:});

