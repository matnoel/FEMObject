function s = calc_epsilon(S,q,varargin)
% function s = calc_epsilon(S,q,varargin)
% S : MODEL
% q : solution

if ~isa(S,'MODEL')
    s = calc_epsilon(q,S,varargin{:});
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
    s = calc_epsilonpc(S,q,varargin{:});
else
    q = unfreevector(S,q);
    if ischarin('node',varargin)
        s = calc_elemfield(S,@epsilonnode,q,varargin{:});
    else
        s = calc_elemfield(S,@epsilon,q,varargin{:});
    end
end

if timeq
    s = TIMEMATRIX(s,T);
end
