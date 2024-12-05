function varargout = pdfplot(upc,varargin)
% function [P,x] = pdfplot(upc,varargin)
% trace la pdf de upc
% appelle [P,x]=plot(upc,'pdf',varargin{:});
%
% See also sampleplot, PCMATRIX/cdfplot, PCMATRIX/plot, PCMATRIX/funplot, RANDVAR/plot
varargout=cell(1,nargout);

if numel(upc)>1
    
    [varargout{:}]=plot(upc,'multipdf',varargin{:});
else
    [varargout{:}]=plot(upc,'pdf',varargin{:});
end    
