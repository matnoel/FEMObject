function Z = getmatrixatstep(at,i)
% function Z = getmatrixatstep(at,i)

s.type = '{}';
s.subs{1} = i;
Z = subsref(at,s);

