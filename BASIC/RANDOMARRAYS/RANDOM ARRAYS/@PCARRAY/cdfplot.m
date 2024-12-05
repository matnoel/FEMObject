function [P,x] = cdfplot(upc,varargin)
% function [P,x] = cdfplot(upc,varargin)
% trace la cdf de upc
% appelle [P,x]=plot(upc,'cdf',varargin{:});

[P,x]=plot(upc,'cdf',varargin{:});