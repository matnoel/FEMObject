function weights = calcStabilisationWeights(assembler)
% weights = stabilisationWeights(assembler)
% Computes weights for SWIP stabilisation operator
% Output is formatted as TSpaceOperators cell array.

Kmax = getConductivityBounds(assembler,'max') ;

weights = calcStabilisationWeights(getOrder(assembler),Kmax,...
    getCellNum(assembler),getTolSVD(assembler)) ;

end