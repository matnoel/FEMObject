function [ur,xi]= random(u,varargin)
% function [ur,xi]= random(u,varargin)

if isa(u.value,'cell')
    if nargin==2
        error('pas prevu')
    end
    ur = u;
    xi = random(RANDVARS(u));
    for i=1:length(u.value)
        if israndom(u.value{i})
            ur.value{i} = randomeval(u.value{i},xi);
        end
    end
    ur = TIMEMATRIX(ur.value,ur.TIMEMODEL,ur.s);
else
    if nargin==1
        ur = u;
        [ur.value,xi] = random(u.value);
        ur = TIMEMATRIX(ur.value,ur.TIMEMODEL,ur.s);
    elseif nargin==2
        ur = u;
        [ur.value,xi] = random(u.value,varargin{1});
        ur = TIMEMATRIX(ur.value,ur.TIMEMODEL,ur.s);
    else
        error('rentrer un ou 2 arguments')
    end
end
