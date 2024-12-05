function disp(u)
% function disp(u)

fprintf(['(RANDPOLY #' num2str(getnumber(u)) ' -> ' class(u) ')'])
disp(' ')
if length(struct2cell(u.param))>0
    disp(u.param)
else
    disp(' ')
end
