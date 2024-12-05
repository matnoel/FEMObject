function varargout = multipdfsampleplot(us,varargin)
% function [P,x] = multipdfsampleplot(us,varargin)
% trace la pdf d'un ensemble de realisations d'un vecteur aleatoire
% appelle [P,x]=multisampleplot(us,'multipdf',varargin{:});
%
% See also sampleplot, cpdfsampleplot

varargout=cell(1,nargout);
[varargout{:}]=sampleplot(us,'multipdf',varargin{:});

