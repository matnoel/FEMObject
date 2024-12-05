function [op,srcOp] = sipOperator(model,K,source,isPeriodic,penalty)
% [op,srcOp] = sipOperator(model,K,source,isPeriodic,penalty)
    
    if nargin<5 || isempty(penalty)
        bounds = mesoBounds(model,K) ;
        penalty = 1.5*penaltyBound(model,bounds,false,false,isPeriodic) ;
    end

    cells = repmat(formatIndex(getOrder(model),getCellNum(model),...
                               (1:getCellNb(model))'),1,2) ;
    diffOp = bilinFormOperator(model,[1 1 0],cells,K) ;
    consOp = consistencyOperator(model,K,isPeriodic) ;
    stabOp = penalty*stabilisationOperator(model,1,isPeriodic) ;

    op = diffOp + stabOp - consOp - consOp' ;

    if nargout > 1
        if ischar(source)
            if strcmp(source,'corrector1')
                coord = getCoord(model) ;
                srcOp = (diffOp-consOp)*coord{1} ;
            elseif strcmp(source,'corrector2')
                coord = getCoord(model) ;
                srcOp = (diffOp-consOp)*coord{2} ;
            end
        else
            srcOp = bilinFormOperator(model,[0 0 0],cells,1)*source ;
        end
    end
end