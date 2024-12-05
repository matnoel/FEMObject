function [i,j,k,val,V] = trimatrixelem(elem,me,s)

numddlelem = getnumddl(elem);
nbddl = getnbddl(elem);
nbelem = getnbelem(elem) ;

numddlelem=reshape(numddlelem',1,1,nbddl,nbelem);

repcol = repmat(permute(numddlelem,[3,2,1,4]),                 [1,nbddl,nbddl,1])+...
         repmat(permute(s(1)*     (numddlelem-1),[1,3,2,4]),   [nbddl,1,nbddl,1])+...
         repmat(        s(1)*s(2)*(numddlelem-1),              [nbddl,nbddl,1,1]);
repcol=  reshape(repcol,nbddl^3,nbelem);

if nargout==4
    % Assemblage des elements
    repu=unique(repcol);
    [temp,repcolu]=ismember(repcol,repu);
    Me=cell(1,nbelem);
    Me(:)={spalloc(length(repu),1,nbddl^3)};


    for e=1:nbelem
    rep=repcolu(:,e);
    Me{e}(rep)=me(:,:,:,e);
    end
    
    [i,j,k]=ind2sub(s,repu);
    val=sum([Me{:}],2);
elseif nargout==5
    % Sans assemblage d'elements
    [i,j,k]=ind2sub(s,repcol(:));
    val=repmat(1:nbelem,[nbddl^3,1]);
    val=val(:);
    V=me(:);
end