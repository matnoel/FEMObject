function varargout = fplot(upc,varargin)
% function [P,x] = fplot(upc,varargin)
% trace la surface de reponse de upc en fonction des variables de base du chaos
% appelle [P,x]=plot(upc,'fun',varargin{:});
%
% See also sampleplot, PCMATRIX/cdfplot, PCMATRIX/pdfplot, RANDVAR/plot

varargout=cell(1,nargout);
[varargout{:}]=plot(upc,'fun',varargin{:});