function S_part = getmodelpart(S,k)
% function S_part = getmodelpart(S,k)
% Get model S_part corresponding to group k of elements of model S
% S_part, S : MODEL

numgroupelem = getnumgroupelemwithparam(S,'partition',k);
S_part = keepgroupelem(S,numgroupelem);
S_part = removenodewithoutelem(S_part);
S_part = keepeleminnode(S_part);

end