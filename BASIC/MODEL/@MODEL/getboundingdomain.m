function D = getboundingdomain(S)
%function D = getboundingdomain(S)
% get the smallest domain that contains the mesh S
% S: MODEL
% D: DOMAIN 

x = getcoord(getnode(S));
m = min(x,[],1)-eps;
M = max(x,[],1)+eps;
D = DOMAIN(size(x,2),m,M);




