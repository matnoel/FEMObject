function P = calc_P_edgeelem(M,Edge1,Edge2)

numNode1 = getnumnode_edgeelem(M,Edge1) ;
numNode2 = getnumnode_edgeelem(M,Edge2) ;

P = calc_P_edges(M,numNode1,numNode2) ;

end

function numNode = getnumnode_edgeelem(M,Edge)

[~,numNodeEdge,~] = intersect(M,Edge) ;
numElem = findelemwithnode(M,numNodeEdge) ;
numNode = [] ;
for p = 1:M.nbgroupelem
    numNode = [numNode ; getnumnode(getelem(M.groupelem{p},numElem,'global'))] ;
end
numNode = unique(numNode) ;

end