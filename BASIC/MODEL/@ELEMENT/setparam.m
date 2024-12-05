function elem = setparam(elem,varargin)
% function elem = setparam(elem,param)
%
% function elem = setparam(elem,paramname,paramval)

if nargin==2
    elem.param = varargin{1};
elseif nargin==3
    elem.param = setfield(elem.param,varargin{1},varargin{2});
end




