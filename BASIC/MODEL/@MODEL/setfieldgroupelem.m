function M=setfieldgroupelem(M,field,value,p)
% function M=setfieldgroupelem(M,field,value,p)
% met le champ field des groupes d'elements p à la valeur value
%
% field et value peuvent etre des cellules pour changer plusieurs champs a
% la fois

if nargin==3
    p = 1:M.nbgroupelem;
end

for i=1:length(p)
    if ~isa(field,'cell')
        M.groupelem{p(i)} = setelementfield(M.groupelem{p(i)},field,value);
    else
        for j=1:length(field)
            M.groupelem{p(i)} = setelementfield(M.groupelem{p(i)},field{j},value{j});
        end
    end
end
