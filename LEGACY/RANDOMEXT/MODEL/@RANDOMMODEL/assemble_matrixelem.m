function M = assemble_matrixelem(S,me,varargin)
% function M = assemble_matrixelem(S,me,varargin)
% s taille de la matrice finale

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
s = [S.nbddl,S.nbddl];
i = [];
j = [];
val = [];
for p=liste
    [ip,jp,valp] = matrixelem(S.groupelem{p},double(me{p}),s);
    i = [i;ip];
    j = [j;jp];
    val = [val;valp];
end

M = sparse(i,j,val,s(1),s(2));


