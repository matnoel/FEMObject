function [op,rhsOp,penalty] = swipOperator(model,K,source,bc,penalty,useSW,useAW)
% [op,rhsOp,penalty] = swipOperator(model,K,source,bc,penalty,useSW,useAW)

if nargin < 7
    useAW = false ;
    if nargin < 6
        useSW = false ;
        if nargin < 5
            penalty = [] ;
            if nargin < 4
                bc = [] ;
                if nargin < 3
                    source = 0 ;
                end
            end
        end
    end
end

% Periodicity and connectivity
periodicity = [false false] ;
if ~isempty(bc)
    periodicity = [isPeriodic(bc,1) isPeriodic(bc,2)] ;
end
connectivity = {mesoConnectivity(model,1,periodicity(1)) ; ...
    mesoConnectivity(model,2,periodicity(2))} ;

% Compute conductivity bounds if necessary
if useAW || useSW || isempty(penalty) || penalty<=0
    bounds = mesoBounds(model,K) ;
end

% Diffusion operator
cells = repmat(formatIndex(getOrder(model),getCellNum(model),...
    (1:getCellNb(model))'),1,2) ;
diffOp = bilinFormOperator(model,[1 1 0],cells,K) ;

% Consistency operator
awOp = [] ;
if useAW
    awOp = weightsOperator(model,bounds(:,4),'average') ;
end
consOp = consistencyOperator(model,K,connectivity,awOp) ;

% Stabilisation operator
if isempty(penalty) || penalty<=0
    penalty = 1.5*penaltyBound(model,bounds,useSW,useAW,periodicity) ;
end
swOp = [] ;
if useSW
    swOp = weightsOperator(model,bounds(:,4),'stabilisation') ;
end
stabOp = penalty*stabilisationOperator(model,1,connectivity,swOp) ;

% BC operators
lBCOp = [] ;
rhsOp = [] ;
if ~isempty(bc)
    missingNitschePenalty = [bc.type]==1 & isempty([bc.factor]) ;
    if any(missingNitschePenalty)
        [bc(missingNitschePenalty).factor] = deal(100*penalty) ; % use same penalty ?
    end
    if nargout > 1
        [lBCOp,rhsOp] = bc.operators(model,K,true) ;
    else
        lBCOp = bc.operators(model,K,true) ;
    end
end

% SWIP operator
op = diffOp + stabOp - consOp - consOp' + lBCOp ;

% Source operator
if nargout > 1
    srcOp = [] ;
    if ischar(source)
        if strcmp(source,'corrector1')
            coord = getCoord(model) ;
            srcOp = -(diffOp-consOp)*coord{1} ;
        elseif strcmp(source,'corrector2')
            coord = getCoord(model) ;
            srcOp = -(diffOp-consOp)*coord{2} ;
        end
    elseif isnumeric(source) || isa(source,'TuckerLikeTensor')
        if norm(source,Inf) < 1e-14 % i.e. source is zero
            % do nothing and leave srcOp empty
        else
            if isnumeric(source)
                if isscalar(source)
                    source = source*TuckerLikeTensor.ones(model.tensorSize);
                else
                    sourceMicro = source ;
                    source = TuckerLikeTensor.ones(model.tensorSize);
                    source.space.spaces{end} = sourceMicro ;
                end
            end
            srcOp = bilinFormOperator(model,[0 0 0],cells,1)*source ;
        end
    else
        error('Unknown source format')
    end
    rhsOp = rhsOp + srcOp ;
    if isempty(rhsOp)
        rhsOp =  TuckerLikeTensor.zeros(model.tensorSize()) ;
    end
end
end