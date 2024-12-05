function disp(PC)
% function disp(PC)

disp(PC.POLYCHAOS,'nopolys')

%if length(PC.groups)>1
for i=1:length(PC.groups)
    fprintf('   group #%d : ',i)
    fprintf('%d ',PC.groups{i});
    fprintf('\n')
end