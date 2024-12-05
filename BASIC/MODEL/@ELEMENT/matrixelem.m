function [i,j,val,V] = matrixelem(elem,me,s,varargin)
% function [i,j,val,V] = matrixelem(elem,me,s,varargin)

numddlelem = getnumddl(elem);
nbddl = getnbddl(elem);

numelem = getnumber(elem);
nbelem = getnbelem(elem) ;

numddlelem = reshape(numddlelem',1,nbddl,nbelem);

repcol = repmat(permute(numddlelem,[2,1,3]),[1,nbddl,1])+...
    repmat(s(1)*(numddlelem-1),[nbddl,1,1]);
repcol = reshape(repcol,nbddl^2,nbelem);


if nargout==3
    % Assemblage des elements
    repu = unique(repcol);
    [temp,repcolu] = ismember(repcol,repu);
    Me = cell(1,nbelem);
    Me(:) = {spalloc(length(repu),1,nbddl^2)};
    
    
    for e=1:nbelem
        rep = repcolu(:,e);
        Me{e}(rep) = me(:,:,e);
    end
    
    [i,j] = ind2sub(s,repu);
    val = sum([Me{:}],2);
elseif nargout==4
    % Sans assemblage d'elements
    [i,j] = ind2sub(s,repcol(:));
    val = repmat(1:nbelem,[nbddl^2,1]);
    val = val(:);
    V = me(:);
end