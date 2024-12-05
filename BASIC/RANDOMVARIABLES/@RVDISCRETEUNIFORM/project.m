function apc = project(a,pc,varargin)
    param = getparam(a);
    apc = PCMATRIX(ones(1,param.Q),[1 1],pc);
end
