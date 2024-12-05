function M=addgroupelem(M,elem)
% function M=addgroupelem(M,elem)

if ~isa(elem,'cell')
elem={elem};
end

for i=1:length(elem)
if ~isa(elem{i},'ELEMENTGEOM')
    error('mauvais type de groupe d''element')
end
M.groupelem{end+1}=elem{i};
end

M.nbgroupelem = length(M.groupelem);
