function varargout = multipdfplot(upc,varargin)
% function [P,x] = multipdfplot(upc,varargin)
% trace la multipdf de upc (doit contenir 2 variables aleatoires)
% appelle [P,x]=plot(upc,'pdf',varargin{:});
%
% See also sampleplot, PCMATRIX/cdfplot, PCMATRIX/plot, PCMATRIX/funplot, RANDVAR/plot
varargout=cell(1,nargout);
[varargout{:}]=plot(upc,'multipdf',varargin{:});
