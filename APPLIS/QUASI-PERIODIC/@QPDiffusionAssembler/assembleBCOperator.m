function [lhs,rhs,time] = assembleBCOperator(assembler)
% [lhs,rhs,time] = assembleBCOperator(assembler)

tic ;

bc = getBC(assembler) ;
microOp = getMicroOperators(assembler) ;
K = getConductivity(assembler) ;
Kmax = getConductivityBounds(assembler) ;
Kmax = Kmax(:,2) ;
[lhs,rhs] = assemble(bc{1},microOp,K,Kmax,1) ;
for i = 2:numel(bc)
    [newLHS,newRHS] = assemble(bc{i},microOp,K,Kmax,1) ;
    lhs = lhs + newLHS ;
    rhs = rhs + newRHS ;
end

% Add constant nullification if only Neumann and periodic BC
NeumannOrPeriodicBC = ismember(cellfun(@getType,bc),[2 5]) ;
if all(NeumannOrPeriodicBC) % only Neumann and Periodic BC
    cstNullOp = assembleCstNullOperator(assembler) ;
    lhs = lhs + cstNullOp ;
end

time = toc ;
end