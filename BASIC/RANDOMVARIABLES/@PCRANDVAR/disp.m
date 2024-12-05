function disp(u)
% function disp(u)

fprintf(['(RANDVAR #' num2str(getnumber(u)) ' -> ' class(u) ')'])
disp(' ')
disp(u.PCMATRIX)