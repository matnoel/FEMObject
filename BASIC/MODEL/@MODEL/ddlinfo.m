function ddlinfo(M,numddl)
%function ddlinfo(M,numddl)
groupddl=getgroupddl(M.node);
if nargin==1
    numddl=(1:getnbddl(M));
end
list = cell(numel(numddl),5);
repblocktemp = ismember(numddl,getddlbloque(M));
repfound = [];
for k=1:numel(numddl)
    for p=1:length(groupddl)
        rep = find(groupddl{p}.numddl==numddl(k));
        [i,j] = ind2sub(size(groupddl{p}.numddl),rep);
        if ~isempty(rep)
            list(k,1) = {num2str(numddl(k))};
            list(k,2) = {num2str(groupddl{p}.node(i))};
            list(k,3) = {getddlname(groupddl{p}.ddlnode,j)};
            list(k,4) = {getddlname(groupddl{p}.ddlnodedual,j)};
            list(k,5) = {num2str(repblocktemp(k))};
            repfound = [repfound,k];
        end
    end
end

fprintf('[numddl,numnode,ddlnode,ddlnodedual,block] = \n')
disp(list(repfound,:))

if numel(numddl)>numel(repfound)
    fprintf('DDL non trouve : \n')
    disp(numddl(repfound));
end
