function weights = calcStabilisationWeights(order,bounds,cellNum,tolSVD)
% weights = calcStabilisationWeights(order,bounds,cellNum,tolSVD)
% Computes SWIP stabilisation operator weights.
% Output is formatted as TSpaceOperators cell array.

if nargin < 4
    tolSVD = 1e-6 ;
    if nargin < 3
        assert(order==2,...
            'Order 3 format requires mesoscopic domain dimensions')
        cellNum = [] ;
    end
end

weights = 2*(bounds*bounds')./bsxfun(@plus,bounds,bounds') ;

expectedNaN = bounds == 0 ;
if any(expectedNaN)
    warning('Some bounds are zero, resulting in NaN weights. Replacing NaN weights by zeros.')
    expectedNaN = find(expectedNaN) ;
    weights(expectedNaN,expectedNaN) = 0 ;
end

if order==2 % then stop here
    weights = {weights} ;
else
    % Order 3 processing (including truncated SVD)
    weights = svdQP(weights,cellNum,tolSVD) ;
end

end