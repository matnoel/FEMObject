function disp(u)
% function disp(u)

disp([' stodim   : ' num2str(u.stodim)])
disp([' storage : ' u.storage])
disp(' ')
if strcmp(u.storage,'PC')
    fprintf(' Polychaos Coefficients ')
    display(u.V)
elseif strcmp(u.storage,'RADIAL')
    fprintf(' Deterministic functions ')
    display(u.V)
    disp(' ')
    
    fprintf(' Stochastic functions ')
    display(u.L)
end


