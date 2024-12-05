function M = splitgroupelem(M,nmax)
% function M = splitgroupelem(M,nmax)
% Separations des groupes d'elements en sous groupes (pour diminuer le temps de calcul)
% nmax : nombre d'element maxi dans le groupe
if nargin==1
    error('definir un nombre maxi d''elements par groupe')
end
for i=1:M.nbgroupelem
elem=M.groupelem{i};
n=getnbelem(elem);
if n>nmax
nb = ceil(n/nmax);    
r=repmat(ceil(n/nb),1,nb);
r(end)=n-sum(r(1:end-1));
r=cumsum(r);
M.groupelem{i}=getelem(elem,1:r(1));
for k=1:length(r)-1
M.groupelem{M.nbgroupelem+k}=getelem(elem,r(k)+1:r(k+1));
end
M.nbgroupelem=M.nbgroupelem+length(r)-1;
end
end

M=changeelemnumber(M);

