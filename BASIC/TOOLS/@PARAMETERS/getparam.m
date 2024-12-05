function param = getparam(P,param,def)
% function param = getparam(P,param,def)

if nargin==1
    param = P.param;
else
    if isfield(P.param,param)
        param = getfield(P.param,param);
    elseif nargin==3
        param = def;
    else
        error('le parametre demande n''existe pas')
    end
end

