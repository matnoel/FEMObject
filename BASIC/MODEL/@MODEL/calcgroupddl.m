function groupddl=calcgroupddl(M)

% CREATION DES GROUPES DE DDL
groupddl=cell(0,1);
for p=1:M.nbgroupelem
groupddlelem = calcgroupddl(M.groupelem{p});
if groupddlelem.nbddl>0
n = length(groupddl)+1;
groupddl{n} = groupddlelem ;
nbddln=groupddl{n}.nbddl;

for j=1:n-1
    nbddlj = groupddl{j}.nbddl ;
    [an,bn] = findddl(groupddl{n}.ddlnode,groupddl{j}.ddlnode)  ;  
    
    if nbddln==nbddlj && length(an)==nbddln
    groupddl{n}.node = union(groupddl{n}.node,groupddl{j}.node)  ;
    groupddl{j}.node=[] ;
    
    elseif  (nbddln<nbddlj && length(an)==nbddln) 
    groupddl{n}.node = setdiff(groupddl{n}.node,groupddl{j}.node);
    
    elseif  (nbddlj<nbddln && length(bn)==nbddlj) 
    groupddl{j}.node = setdiff(groupddl{j}.node,groupddl{n}.node);
    
    else
    error('non programme')
    end
 
end
        
end
end


repelim = [];
for p=1:length(groupddl)
if isempty(groupddl{p}.node)
repelim=[repelim,p];    
end
end

groupddl(repelim)=[];

