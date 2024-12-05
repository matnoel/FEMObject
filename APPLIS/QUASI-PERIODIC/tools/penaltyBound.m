function penalty = penaltyBound(bounds,traceCst,cellSize,useSW,useAW,cellConn)
% penalty = penaltyBound(bounds,traceCst,cellSize,useSW,useAW,cellConn)

if nargin < 6
    if nargin < 5
        useAW = false ;
        if nargin < 4
            useSW = false ;
        end
    end 
    assert(~(useAW||useSW),['Cannot apply required weights: ',...
        'mesoscopic connectivity missing.']) ;
end

% Weights (if SWIP)
if useAW
    maxAW = calcAverageWeights(2,bounds(:,4)) ; % force order 2 format
    maxAW = maxAW{1}.*cellConn ;
    maxAW = full(max(maxAW(maxAW>0))) ;
    if isempty(maxAW); maxAW=1 ; end % robustness (case 1 cell)
else
    maxAW = .5 ;
end
if useSW
    minSW = calcStabilisationWeights(2,bounds(:,4)) ;
    minSW = minSW{1}.*cellConn ;
    minSW = full(min(minSW(minSW>0))) ;
    if isempty(minSW); minSW=1 ; end % robustness (case 1 cell)
else
    minSW = 1;
end

% Various values
faceMax = max(cellSize) ;
kMinusMin = min(bounds(:,1)) ;
kPlusMax = max(bounds(:,4)) ;

% Compute penalty lower bound and apply factor
penalty = 4*faceMax*(kPlusMax*traceCst*maxAW)^2/(kMinusMin*minSW) ;
end