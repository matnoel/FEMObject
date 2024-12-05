function ur = randomeval(u,varargin)
% function ur = randomeval(u,varargin)

ur = u;
if isa(u.value,'cell')
    for i=1:length(u.value)
        if israndom(u.value{i})
            ur.value{i} = randomeval(u.value{i},varargin{:});
        end
        ur = TIMEMATRIX(ur.value,ur.TIMEMODEL,ur.s);
    end
else
    ur.value = randomeval(u.value,varargin{:});
    ur = TIMEMATRIX(ur.value,ur.TIMEMODEL,ur.s);
end
