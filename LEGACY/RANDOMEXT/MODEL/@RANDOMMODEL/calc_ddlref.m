function rep = calc_ddlref(M,Mrep,ddl)
% function rep = calc_ddlref(M,Mrep,ddl)

xnode = getcoord(M.node);
DX = zeros(M.nbnode,0);
for j=1:length(Mrep)
    DX = [DX,xnode(:,j)-Mrep(j)];
end

normDX = sqrt(sum(DX.^2,2));
numrep = find(normDX==min(normDX));
numnode = numrep(1);

numrep = findddl(M.groupddl{M.repnodeingroupddl(numnode,1)}.ddlnode,ddl);
rep = M.groupddl{M.repnodeingroupddl(numnode,1)}.numddl(M.repnodeingroupddl(numnode,2),numrep);
