function P = setparam(P,param,value,varargin)
% function P = setparam(P,param,value,varargin)

if nargin>=3
    P.param = setfield(P.param,param,value);
    
    if length(varargin)>0
        P =  setparam(P,varargin{:});
    end
end


