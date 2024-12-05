function varargout = cdfplot(u,varargin)
% function [P,x] = pdfplot(u,varargin)
% trace la pdf de u
% appelle [P,x]=plot(u,'pdf',varargin{:});
varargout=cell(1,nargout);
[varargout{:}]=plot(u,'cdfpdf',varargin{:});