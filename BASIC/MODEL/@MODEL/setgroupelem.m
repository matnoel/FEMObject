function M=setgroupelem(M,p,elem)
% function M=setgroupelem(M,p,elem)

if nargin==2
M.groupelem = p;
else
M.groupelem{p}=elem;
end
M.nbgroupelem = length(M.groupelem);
