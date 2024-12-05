function S = divideelem(M)
% function S = divideelem(M)

S = MODEL(M.mode);
S = addnode(S,M.node);

for i=1:M.nbgroupelem
    elem = M.groupelem{i};
    [subelem,nodeplus] = divideelem(elem,S.node);
    S = addnode(S,nodeplus);
    S = addelem(S,subelem);
end

S = unique(S);
S = changenodenumber(S,[1:S.nbnode]);
S = changeelemnumber(S);