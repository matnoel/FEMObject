function varargout = cdfplot(upc,varargin)
% function [P,x] = cdfplot(upc,varargin)
% trace la cdf de upc
% appelle [P,x]=plot(upc,'cdf',varargin{:});
%
% See also sampleplot, PCMATRIX/pdfplot, PCMATRIX/funplot, RANDVAR/plot

varargout=cell(1,nargout);
[varargout{:}]=plot(expand(upc),'cdf',varargin{:});