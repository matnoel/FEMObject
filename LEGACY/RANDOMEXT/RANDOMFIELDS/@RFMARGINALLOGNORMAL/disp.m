function disp(ma)
% function disp(ma)

disp('  Parameters of the marginal distribution')
if ~ma.selfstat
    disp('  (mu and sigma are the mean and standard deviation of the associated gaussian distribution)')
end
a = struct('mu',ma.mu,'sigma',ma.sigma,'x0',ma.x0);
disp(a);
