function varargout = pdfsampleplot(us,varargin)
% function [P,x] = pdfsampleplot(us,varargin)
% trace la pdf d'un ensemble de realisations
% appelle [P,x]=sampleplot(us,'pdf',varargin{:});
%
% See also sampleplot, cpdfsampleplot

varargout=cell(1,nargout);
if any(size(us)==1)
    [varargout{:}]=sampleplot(us,'pdf',varargin{:});
else
    [varargout{:}]=sampleplot(us,'multipdf',varargin{:});    
end

