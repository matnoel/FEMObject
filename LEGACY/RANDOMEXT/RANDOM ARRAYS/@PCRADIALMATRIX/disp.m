function disp(u)
% function disp(u)

disp([' -> number of radial functions : m = ' num2str(u.m)])
disp([' -> stochastic dimension : P = ' num2str(getP(u.L))])
if iscell(u.V)
    disp(' -> storage cell')
else
    disp(' -> storage double')
end
disp(' ')
% disp(' -> stochastic functions : ')
% disp(u.L);
