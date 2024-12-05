function T = setapproxparam(T,r,e,varargin)
% function p = setapproxparam(T,r,e)

T.approxparam = setparam(T.approxparam,r,e);
if length(varargin)>=2
    T = setapproxparam(T,varargin{:});
end

