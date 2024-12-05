function A = subsref(a,s)
% function A = subsref(a,s)

A = subsref(a.forms{1},s);
for i=2:length(a.forms)
    A = A + subsref(a.forms{i},s);
end


