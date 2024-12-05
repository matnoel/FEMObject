function elem = addparam(elem,varargin)
% function elem = addparam(elem,param)
%
% function elem = addparam(elem,paramname,paramval)

if nargin==2
    error('pas programme')
    elem.param = varargin{1};
elseif nargin==3
    elem.param=setfield(elem.param,varargin{1},varargin{2});
end




