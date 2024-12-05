function disp(u)
% function disp(u)

fprintf('(RANDVARFUNCTION)')
disp(' ')
if isa(u.fun,'inline')
    arg = symvar(u.fun);
    disp(u.fun)
else
    arg = symvar(func2str(u.fun));
    disp(' ')
    disp('    Function_handle :')
    disp(u.fun)
end
for i=1:length(u.param)
    fprintf([arg{i} ' = '])
    disp(u.param{i})
end
