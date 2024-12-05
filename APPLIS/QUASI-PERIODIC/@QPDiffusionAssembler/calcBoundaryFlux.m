function flux = calcBoundaryFlux(assembler,v,bCells)
% flux = calcBoundaryFlux(assembler,v,bCells)
% Shortcut to call on QPModel homonym method.

flux = calcBoundaryFlux(getModel(assembler),v,...
    getMicroOperators(assembler),bCells,getConductivity(assembler)) ;
end