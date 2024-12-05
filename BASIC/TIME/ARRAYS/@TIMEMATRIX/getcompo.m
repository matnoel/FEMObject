function Z = getcompo(at,i)
% function Z = getcompo(at,i)

s.type = '()';
s.subs{1} = i;
Z = subsref(at,s);

