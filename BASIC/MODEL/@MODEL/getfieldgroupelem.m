function value=getfieldgroupelem(M,field,p)
% function value=getfieldgroupelem(M,field,p)
% donne la valeur du champ field du groupe d'elements p
% si p est un vecteur, value est une cellule
%

if nargin==2
    p = 1:M.nbgroupelem;
end

value = cell(1,length(p));
for i=1:length(p)
    value{i} = getelementfield(M.groupelem{p(i)},field);
end

if length(value)==1
    value=value{1};
end


