function disp(u)
% function disp(u)

fprintf(['(MATERIAL #' num2str(getnumber(u)) ' -> ' class(u) ')'])
disp(' ')
for k=1:size(u.param,1)
    if isa(u.param{k,2},'RANDVAR')
        u.param{k,2} = [class(u.param{k,2}) ' #' num2str(getnumber(u.param{k,2}))];
    end
end
disp(u.param)