function [c2pElem,c2pRefNodes,c2pRefEdge] = find_edges(connecChild,connecParent)
% function [c2pElem,c2pRefNodes,c2pRefEdge] = find_edges(connecChild,connecParent)
% Locate edges from connecChild in connecParent

nbElemParent = size(connecParent,1) ;
nbNodeElemParent = size(connecParent,2) ;

EdgesList = [] ;
for i=1:(nbNodeElemParent-1)
    EdgesList = [EdgesList ; [connecParent(:,i) connecParent(:,i+1)]] ;
end
EdgesList = [EdgesList ; [connecParent(:,end) connecParent(:,1)]] ;

[~,Loc] = ismember(connecChild,EdgesList,'rows');
if any(Loc==0) % some parents not found, try to swap 1st and 2nd nodes
    [~,Loc2] = ismember(connecChild,circshift(EdgesList,[0 1]),'rows');
    Loc = Loc + Loc2 ; % non-zero values are on different lines
end

c2pElem = 1+mod(Loc-1,nbElemParent) ;
c2pRefEdge = 1+fix((Loc-1)/nbElemParent) ;
c2pRefNodes = [c2pRefEdge 1+mod(c2pRefEdge,nbNodeElemParent)] ;
if exist('Loc2','var') && ~isempty(Loc2) % Swap for those from Loc2
    c2pRefNodes(Loc2~=0,:) = [c2pRefNodes(Loc2~=0,2) c2pRefNodes(Loc2~=0,1)];
end
    
% restore zeros
Loc_zeros = find(Loc==0) ;
c2pElem(Loc_zeros) = 0 ;
c2pRefEdge(Loc_zeros) = 0 ;
c2pRefNodes(Loc_zeros,:) = 0 ;

end