function v = getsubmatrix(v,ind,dim)
% function u = getsubmatrix(v,ind,dim)
% ind = indices a garder suivant la dimension dim
% Par exemple : si v a 4 dimensions et que dim=2
% u=v(:,ind,:,:)
% si n=length(dim)>1, ind doit avoir n colonnes
% Par exemple : si v a 4 dimensions et que dim=[2,3]
% u=v(:,ind(:,1),ind(:,2),:)

if length(dim)==1
    ind=ind(:);
end

rep = cell(1,ndims(v));
rep(:) = {':'};
for i=1:length(dim)
    rep{dim(i)} = ind(:,i);
end

v.double = v.double(rep{:});
