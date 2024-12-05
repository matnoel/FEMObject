function A=fundim(A,fun,dim)
% Application de la fonction 'fun' sur l'ensemble des
% dimensions dim.
% ATTENTION : en HSEP, dim correspond a la dimensions absolue

% Reperer la dim dans la 
for i=1:length(dim)
    recalage=0;
    for k=1:getvar2dim(A.tree,dim(i),1)-1
    recalage=recalage+totaldim(A.F{1,k});
    end
    A.F(:,getvar2dim(A.tree,dim(i),1))=cellfun(@(a) fundim(a,fun,dim(i)-recalage) ,...
                       A.F(:,getvar2dim(A.tree,dim(i),1)),'UniformOutput',0);
end