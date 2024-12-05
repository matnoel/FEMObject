function [errorVal,diff] = compare2FE(model,tensQP,matFE,feCoord)
% [errorVal,diff] = compareResidual2FE(model,tensQP,matFE,feCoord)

if nargin>3 && ~isempty(feCoord)
    qpCoord = smooth(model,getDomainCoord(model)) ;
    matFE = relativeSort(matFE,feCoord,qpCoord) ;
end

tensFE = tensorize(model,matFE) ;
diff = tensQP-tensFE ;
errorVal = norm(diff)/norm(tensFE) ;

end