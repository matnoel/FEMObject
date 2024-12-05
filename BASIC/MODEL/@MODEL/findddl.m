function [listddl,ddlsto] = findddl(M,ddls,varargin)
% function listddl = findddl(M,ddls,numnode)
% function listddl = findddl(M,ddls,D)
% D : object de type POINT ou NODE
% ddls : Objet DDL ou cell ou char
% ddls = 'all' cherche tous les ddls
%
% function listddl = findddl(M,ddls,D,'free')
% donne le reperage de ce ddl dans l'ensemble des ddl libres

[rep1,pos1] = isclassin('GEOMOBJECT',varargin);
[rep2,pos2] = isclassin('NODE',varargin);
[rep3,pos3] = isclassin('double',varargin);
if rep1
    D = varargin{pos1};
    numnode = ispointin(D,POINT(M.node));
    numnode = getnumber(M.node,numnode);
elseif rep2
    D = POINT(varargin{pos2});
    numnode = intersect(M,D);
    numnode = getnumber(M.node,numnode);
elseif rep3
    numnode = varargin{pos3};
else
    numnode = getnumber(M.node);
end


listddl = zeros(length(ddls),length(numnode));
groupddl = getgroupddl(M.node);
ddlsto = DDL();
for j=1:length(groupddl)
    if ischarin('dual',varargin)
        [rep1,rep2] = findddl(groupddl{j}.ddlnodedual,ddls);
        ddltemp = getddl(groupddl{j}.ddlnodedual);
    else
        [rep1,rep2] = findddl(groupddl{j}.ddlnode,ddls);
        ddltemp = getddl(groupddl{j}.ddlnode);
    end
    [num,repnode1,repnode2] = intersect(groupddl{j}.node,numnode);
    listddl(rep2,repnode2) = groupddl{j}.numddl(repnode1,rep1)';
    if ~isempty(repnode1) && getnbddl(ddltemp)==0
        ddlsto = addddl(ddltemp);
    end
end

listddl = nonzeros(listddl(:));

if ischarin('free',varargin)
    ddlfree = getddlfree(M.BCOND);
    [a,listddl] = ismember(listddl,ddlfree(:));
    listddl = listddl(find(listddl));
end
