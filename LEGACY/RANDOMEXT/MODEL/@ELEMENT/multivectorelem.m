function [i,val,n,nmat] = multivectorelem(elem,fe)

numddlelem = getnumddl(elem);
numelem = getnumber(elem);
nbelem = getnbelem(elem) ;

n=size(fe,2);
nmat=size(fe,4);
fe = permute(fe,[1,2,4,3]);
fe = reshape(fe,[size(fe,1),size(fe,2)*size(fe,3),size(fe,4)]);
replig=numddlelem;
repu = unique(replig);
[temp,repligu]=ismember(replig,repu);

Fe=cell(1,nbelem);
Fe(:)={sparse(length(repu),n*nmat)};
for e=1:nbelem
rep=repligu(e,:);
Fe{e}(rep,:)=fe(:,:,e);
end
val = reshape([Fe{:}],length(repu)*n*nmat,nbelem);
val = reshape(sum(val,2),length(repu),n*nmat);
i=repu;
