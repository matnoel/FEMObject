function wOp = weightsOperator(model,bounds,type)
% wOp = weightsOperator(model,bounds,type)
% Bounds must be upper and as from mesoBounds.
% Type is either 1 or 'stabilisation' for stabilisation weights,
% or 2 or 'average' for average weights.

    if isempty(bounds)
        bounds = ones(getCellNb(model),1) ;
    end

    % Calculate weights
    order = getOrder(model) ;
    switch type
      case {1,'stabilisation'}
        w = calcStabilisationWeights(order,bounds,getCellNum(model), ...
                                      getTolSVD(model));
      case {2,'average'}
        w = calcAverageWeights(order,bounds,getCellNum(model), ...
                                      getTolSVD(model));
      otherwise
        error('Unknown weight type.')
    end


    % Space
    wSpace = [w ; {ones(getNbCellDoF(model))}] ;
    wSpace = TSpaceOperators(wSpace) ;

    % Core
    if order == 2
        wCore = DiagonalTensor(1,2) ;
    else
        wCore = FullTensor(eye(wSpace.dim(1)),order,wSpace.dim) ;
    end

    % Tensor
    wOp = TuckerLikeTensor(wCore,wSpace) ;

end