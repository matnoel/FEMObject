function disp(ma)
% function disp(ma)

disp(['(' class(ma) ')'])
disp(' ')

disp([' randvar type : '  class(ma.RV)])
if length(fieldnames(ma.param))>0
    disp([' with parameters'])
    disp(ma.param)
else
    disp(' ')
end
