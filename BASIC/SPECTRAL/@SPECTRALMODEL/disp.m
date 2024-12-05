function disp(S)
% function disp(S)

fprintf('dim = %d\n',S.dim)
fprintf('number of points in dimension 1 = %d\n',getnbpoints(S.L));
if S.dim==2
    fprintf('number of points in dimension 2 = %d\n',getnbpoints(S.L));
end
fprintf('Boundary conditions \n')
display(S.bc)

