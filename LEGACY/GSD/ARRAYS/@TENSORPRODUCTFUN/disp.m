function disp(u)
% function disp(u)

disp(['number of dimensions : m = ' num2str(u.m)])
disp(['dimensions : ' num2str(getdim(u.PRODUCTSPACE))])
disp(['factor : '])
disp(u.factor)
disp(['functions : '])
display(u.phi);

