function S = separateelemwithnbddl(S)
% function S =  separateelemwithnbddl(S)
% fait en sorte que chaque tous les elements d'un groupe 
% d'elements aient le meme nombre de ddl

node = getnode(S);
nbddlpernode = getnbddlpernode(node);

if length(unique(nbddlpernode))>1

for p=1:getnbgroupelem(S)
   elem = getgroupelem(S,p);
 
   if isenrich(elem) || strcmp(getoption(elem),'FACE')
   connec = calc_conneclocal(elem,node);
   nbddlpernode = getnbddlpernode(node);
   connecnbddl = reshape(nbddlpernode(connec),size(connec));
   
   nbddlperelem = sum(connecnbddl,2);
   N = unique(nbddlperelem);  
   if length(N)>1
       newelem = cell(1,length(N));
   for k=1:length(N)
      rep = find(nbddlperelem==N(k));
      newelem{k} = getelem(elem,rep);
   end
   
   S = setgroupelem(S,p,newelem{1});
   S = addgroupelem(S,newelem(2:end));
   end
   end

end
end

