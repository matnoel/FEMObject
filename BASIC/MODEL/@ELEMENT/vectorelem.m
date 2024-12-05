function [i,val,n] = vectorelem(elem,fe)
% function [i,val,n] = vectorelem(elem,fe)

numddlelem = getnumddl(elem);
% numelem = getnumber(elem);
nbelem = getnbelem(elem) ;
nbddl = getnbddl(elem);
n = size(fe,2);

replig = numddlelem;
repu = unique(replig);
[temp,repligu] = ismember(replig,repu);
Fe = cell(1,nbelem);
Fe(:) = {spalloc(length(repu),n,nbddl*n)};

repligun = repmat(repligu(:),1,n) + ones(numel(repligu),1)*length(repu)*[0:n-1];
repligun = reshape(repligun,size(repligu,1),nbddl*n);


for e=1:nbelem
    rep = repligun(e,:);
    Fe{e}(rep) = fe(:,:,e);
end


val = reshape([Fe{:}],length(repu)*n,nbelem);
val = reshape(sum(val,2),length(repu),n);
i = repu(:);
val = sparse(val);

