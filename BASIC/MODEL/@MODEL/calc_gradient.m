function s = calc_gradient(S,q,varargin)
% function s = calc_gradient(S,q)
% S: MODEL
% q: solution

if israndom(S) || israndom(q)
    s = calc_epsilonpc(S,q,varargin{:});
else 
    q = unfreevector(S,q);
    if ischarin('node',varargin)
        s = calc_elemfield(S,@gradientnode,q,varargin{:});
    else
        s = calc_elemfield(S,@gradient,q,varargin{:});
    end
end
