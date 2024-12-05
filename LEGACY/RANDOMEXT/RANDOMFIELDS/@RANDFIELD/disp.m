function disp(u)
% function disp(u)

disp(['(' class(u) ')'])

if ~u.selfcorrel
    disp('  Nonlinear functional of a normalized gaussian random field')
end

disp([' - MARGINAL DISTRIBUTION : (' class(u.marginal) ')'])
disp(u.marginal)

fprintf(' - CORRELATION ')
if ~u.selfcorrel
    fprintf('(associated with the gaussian random field) ')
end
disp(u.correl)

