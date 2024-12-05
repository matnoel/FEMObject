function disp(u)
% function disp(u)

fprintf(['(RANDVAR #' num2str(getnumber(u)) ' -> ' class(u) ')'])
disp(' ')
param = u.param;
for i=1:size(param,1)
    if israndom(param{i,2})
        param{i,2} = num2str(param{i,2});
    end
end
param = param';

disp(struct(param{:}))