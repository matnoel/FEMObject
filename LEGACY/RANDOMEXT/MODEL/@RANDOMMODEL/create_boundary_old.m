function [frontiere,interfaceglob,interface] = create_boundary(M)

if M.nbgroupelem==1
[elembord,nodebord] = boundary(M.groupelem{1},M.node);
elembord = copy_ddl(elembord,M.groupelem{1});
frontiere=MODEL(M.mode);
frontiere=addnode(frontiere,nodebord);
frontiere=addelem(frontiere,elembord);

else
[elembord,nodebord] = boundary(M.groupelem{1},M.node);
elembord = copy_ddl(elembord,M.groupelem{1});
    bord{1}=MODEL(M.mode);
    bord{1}=addnode(bord{1},nodebord);
    bord{1}=addelem(bord{1},elembord);

bordelem = bord{1};

for i=2:M.nbgroupelem
[elembord,nodebord]= boundary(M.groupelem{i},M.node);
elembord = copy_ddl(elembord,M.groupelem{i});
    bord{i}=MODEL(M.mode);
    bord{i}=addnode(bord{i},nodebord);
    bord{i}=addelem(bord{i},elembord);
    bordelem=union(bordelem,bord{i},'norenum') ;
end


for i=1:length(bord)
for j=i+1:length(bord)
interface{i}{j} = intersect(bord{i},bord{j});

if ~exist('interfaceglob')
interfaceglob = interface{i}{j} ;
else
interfaceglob = union(interfaceglob,interface{i}{j},'norenum');
end    

end
end

frontiere = setdiff(bordelem,interfaceglob);
end

frontiere = createddlnode(frontiere,M.node);
frontiere.BCOND = M.BCOND;
