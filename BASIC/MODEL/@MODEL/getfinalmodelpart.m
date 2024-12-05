function S_part = getfinalmodelpart(S,k)
% function S_part = getfinalmodelpart(S,k)
% Get model S_part corresponding to group k of elements of model S
% Finalize model S_part
% Transfer boundary conditions from model S to model S_part
% S_part, S : MODEL

S_part = getmodelpart(S,k);
% S_part = createddlnode(S_part);
S_part = final(S_part);
S_part = transfercl(S,S_part);

end