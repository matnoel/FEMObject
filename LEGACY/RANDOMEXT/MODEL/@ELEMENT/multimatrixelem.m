function [i,j,val,nmat] = multimatrixelem(elem,me,S)

numddlelem = getnumddl(elem);
nbddl = getnbddl(elem);

numelem = getnumber(elem);
nbelem = getnbelem(elem) ;
numddlelem=reshape(numddlelem',1,nbddl,nbelem);

repcol = repmat(permute(numddlelem,[2,1,3]),[1,nbddl,1])+...
             repmat(S.nbddl*(numddlelem-1),[nbddl,1,1]);
repcol=reshape(repcol,nbddl^2,nbelem);
repu=unique(repcol);
[temp,repcolu]=ismember(repcol,repu);

nmat=size(me,4);

if nmat==1
Me=cell(1,nbelem);
Me(:)={spalloc(length(repu),1,elem.nbddl^2)};
for e=1:nbelem
rep=repcolu(:,e);
Me{e}(rep)=me(:,:,e);
end 
s = [S.nbddl,S.nbddl];
[i,j]=ind2sub(s,repu);
val=sum([Me{:}],2);
else
    
i=[];
j=[];
val=[];

for k=1:nmat
Me=cell(1,nbelem);
Me(:)={spalloc(length(repu),1,elem.nbddl^2)};
for e=1:nbelem
rep=repcolu(:,e);
Me{e}(rep)=me(:,:,e,k);
end 
i = [i;repu(:)];
j = [j;repmat(k,length(repu),1)];
val=[val;sum([Me{:}],2)];
end

end