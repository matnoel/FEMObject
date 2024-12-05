function varargout = pdfplot(upc,varargin)
% function [P,x] = pdfplot(upc,varargin)
% trace la pdf de upc
% appelle [P,x]=plot(upc,'pdf',varargin{:});

[P,x]=plot(upc,'pdf',varargin{:});
