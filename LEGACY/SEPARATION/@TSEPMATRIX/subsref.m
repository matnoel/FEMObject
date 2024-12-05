function v = subsref(A,s)
% function v = subsref(A,s)

if strcmp(s(1).type,'.')
    switch s(1).subs
        case 'alpha'
            if length(s)==1
                v = A.alpha;
            else
                v = A.alpha(s(2).subs{:});
            end
        case 'F'
            if length(s)==1
                v = A.F;
            else
                v = subsref(A.F,s(2:end));
            end
        case 'dim'
            v = A.dim;
        case 'm'
            if length(s)==1
                v = A.m;
            else
                v = A.m(s(2).subs{:});
            end
    end
else
    error('not implemented')
end
