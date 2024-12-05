function [i,j,val,V] = mixedmatrixelem(elem1,elem2,me,s,varargin)

numddlelem1 = getnumddl(elem1);
nbddl1 = getnbddl(elem1);
numddlelem2 = getnumddl(elem2);
nbddl2 = getnbddl(elem2);

numelem1 = getnumber(elem1);
nbelem1 = getnbelem(elem1) ;
numelem2 = getnumber(elem2);
nbelem2 = getnbelem(elem2) ;

numddlelem1=reshape(numddlelem1',1,nbddl1,nbelem1);
numddlelem2=reshape(numddlelem2',1,nbddl2,nbelem2);

repcol = repmat(permute(numddlelem1,[2,1,3]),[1,nbddl2,1])+...
             repmat(s(1)*(numddlelem2-1),[nbddl1,1,1]);
repcol=reshape(repcol,nbddl1*nbddl2,nbelem1);


if nargout==3
    % Assemblage des elements
    repu=unique(repcol);
    [temp,repcolu]=ismember(repcol,repu);
    Me=cell(1,nbelem1);
    Me(:)={spalloc(length(repu),1,nbddl1*nbddl2)};


    for e=1:nbelem1
    rep=repcolu(:,e);
    Me{e}(rep)=me(:,:,e);
    end 

    [i,j]=ind2sub(s,repu);
    val=sum([Me{:}],2);
elseif nargout==4
    % Sans assemblage d'elements
    [i,j]=ind2sub(s,repcol(:));
    val=repmat(1:nbelem1,[nbddl1*nbddl2,1]);
    val=val(:);
    V=me(:);
end