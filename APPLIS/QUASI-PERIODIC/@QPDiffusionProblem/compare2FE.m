function [errorVal,diff] = compare2FE(pb,tensQP,matFE,errorIndicator,varargin)
% [errorVal,diff] = compare2FE(pb,tensQP,matFE,errorIndicator,varargin)

model = getModel(pb) ;

if ischarin('feModel',varargin)
    feModel = getcharin('feModel',varargin) ;
    feCoord = getcoord(getnode(feModel)) ;
    if length(feCoord)>length(matFE) % fixed DoF are missing
        matFE = unfreevector(feModel,matFE) ;
    end
    qpCoord = smooth(model,getDomainCoord(model)) ;
    matFE = relativeSort(matFE,feCoord,qpCoord) ;
end

matQP = untensorize(model,tensQP) ;
diff = matQP - matFE ;

if ischar(errorIndicator) && strcmp(errorIndicator,'residual')
    lhs = getLHSOperator(pb) ;
    tensFE = tensorize(model,matFE) ;
    diffT = tensQP-tensFE ;
    compressor = Truncator('tolerance',getTolSVD(pb)) ;
    diffT = truncate(compressor,diffT) ;
    lDiffT = truncate(compressor,lhs*diffT) ;
    errorVal = sqrt(abs(dot(lDiffT,diffT)))/norm(tensFE) ;
elseif isa(errorIndicator,'function_handle')
    errorVal = errorIndicator(matFE,matQP) ;
else
    error('errorIndicator not recognized.') ;
end

end