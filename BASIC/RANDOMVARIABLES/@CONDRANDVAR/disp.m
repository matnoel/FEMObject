function disp(u)
% function disp(u)

paramnames = cell(length(u.funparam),1);
paramnames(:) = {'a'};
for i=1:length(u.funparam)
    paramnames{i} = strcat(paramnames{i},num2str(i));
end
paramarg = ['(' paramnames{1}];
for i=2:length(paramnames)
    paramarg = [paramarg ',' paramnames{i}];
end
paramarg = [paramarg ')'];

fprintf(['(CONDRANDVAR #' num2str(getnumber(u)) ' -> ' func2str(u.X) paramarg ')'])
disp(' ')
param = u.funparam;
param = cell(2,length(u.funparam));
param(2,:) = u.funparam;
param(1,:) = paramnames;

disp(struct(param{:}))

for i=1:length(u.funparam)
    if isa(u.funparam{i},'inline') || isa(u.funparam{i},'function_handle')
        arg = symvar(func2str(u.funparam{i}));
    end
end
if exist('arg')
    for i=1:length(u.Y)
        fprintf('%s = ',arg{i})
        disp(u.Y{i});
    end
end

