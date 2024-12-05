function s = calc_sigma(S,q,varargin)
% function s = calc_sigma(S,q)
% S: MODEL
% q: solution

if ~isa(S,'MODEL')
    s = calc_sigma(q,S,varargin{:});
    return
end

if isa(q,'TIMEMATRIX')
    timeq = 1;
    T = gettimemodel(q);
    q = getvalue(q);
else
    timeq = 0;
end

if israndom(S) || israndom(q)
    s = calc_sigmapc(S,q,varargin{:});
else
    q = unfreevector(S,q);
    if ischarin('node',varargin)
        s = calc_elemfield(S,@sigmanode,q,varargin{:});
    else
        if isa(S,'LSMODEL') && getnblevelsets(S)>0
            s = calc_elemfield(S,@sigmals,getlevelsets(S),q,varargin{:});
        else
            s = calc_elemfield(S,@sigma,q,varargin{:});
        end
    end
end

if timeq
    s = TIMEMATRIX(s,T);
end
