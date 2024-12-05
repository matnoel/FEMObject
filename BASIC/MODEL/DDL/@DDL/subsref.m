function v = subsref(u,s)
% function v = subsref(u,s)

v = getddl(u,s.subs{1});