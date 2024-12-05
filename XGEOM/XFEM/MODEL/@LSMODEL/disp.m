function disp(S)
% function disp(S)

disp(S.MODEL)

if length(S.ls)>0
    fprintf('MODEL''S (LEVELSETS)\n')
    disp(S.ls,'novalue')
end

