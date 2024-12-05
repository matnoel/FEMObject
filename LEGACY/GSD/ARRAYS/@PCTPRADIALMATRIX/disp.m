function disp(u)
% function disp(u)

disp([' -> number of radial functions : m = ' num2str(u.m)])
if iscell(u.V)
    disp(' -> storage cell')
else
    disp(' -> storage double')
end
disp(' ')
