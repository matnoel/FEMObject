function weights = calcAverageWeights(assembler)
% weights = calcAverageWeights(assembler)
% Computes weights for SWIP averages.
% Output is formatted as TSpaceOperators cell array.

Kmax = getConductivityBounds(assembler,'max') ;

weights = calcAverageWeights(getOrder(assembler),Kmax,...
    getCellNum(assembler),getTolSVD(assembler)) ;

end