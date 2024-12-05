function varargout = pdfplot(u,varargin)
% function [P,x] = pdfplot(u,varargin)
% trace la pdf de u
% appelle [P,x]=plot(upc,'pdf',varargin{:});
%
% See also sampleplot, cdfplot, plot, funplot, plot
varargout=cell(1,nargout);

if israndom(u) && numel(u)>1   
    [varargout{:}]=statplot(u,'multipdf',varargin{:});
else
    [varargout{:}]=statplot(u,'pdf',varargin{:});
end    
