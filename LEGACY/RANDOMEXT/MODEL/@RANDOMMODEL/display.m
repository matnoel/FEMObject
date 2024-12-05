function display(S)
% function display(S)

disp(' ')
disp([inputname(1) ' = (' class(S) ')'])
disp(' ')
v = struct('mode',S.mode,'nbelem',S.nbelem,'nbnode',S.nbnode,'nbgroupelem',S.nbgroupelem,'nbddl',S.nbddl);
disp(v)

if israndom(S)
    fprintf('MODEL''S (RANDVARS) \n')
    disp(RANDVARS(S))
end
if length(S.ls)>0
    fprintf('MODEL''S (LEVELSETS)\n')
    disp(S.ls,'novalue')
end

fprintf('MODEL''S (MATERIALS)\n')
disp(MATERIALS(S))

