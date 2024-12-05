function T = setevolparam(T,r,e,varargin)
% function p = setevolparam(T,r,e)

T.evolparam = setparam(T.evolparam,r,e);
if length(varargin)>=2
    T = setevolparam(T,varargin{:});
end

