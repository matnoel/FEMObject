function penalty = calcPenalty(assembler)

factor = 1.01 ; % Arbitrary provided factor > 1

% Conductivity bounds
Kbounds = getConductivityBounds(assembler) ;
if isempty(Kbounds)
    cA = getConductivityAssembler(assembler) ;
    Kbounds = calcConductivityBounds(cA) ;
    cA = setConductivityBounds(cA,Kbounds) ;
    assembler = setConductivityAssembler(assembler,cA) ;
end
kMinusMin = min(Kbounds(:,1)) ;
kPlusMax = max(Kbounds(:,2)) ;

% Weights (if SWIP)
useAW = getUseAverageWeights(assembler) ;
useSW = getUseStabilisationWeights(assembler) ;
if useSW || useAW
    % Get bidirectional cell connectivity
    cellConn = mesoConnectivity(2,getCellNum(assembler),1, ...
                                       isPeriodic(assembler,1)) ;
    cellConn = cellConn | mesoConnectivity(2,getCellNum(assembler),2, ...
                                           isPeriodic(assembler,2)) ;
end
if useAW
    aWmax = calcAverageWeights(2,Kbounds(:,2)) ; % force order 2 format
    aWmax = aWmax{1}.*cellConn ;
    aWmax = full(max(aWmax(aWmax>0))) ;
    if isempty(aWmax); aWmax=1 ; end % robustness (case 1 cell)
else
    aWmax = .5 ;
end
if useSW
    sWmin = calcStabilisationWeights(2,Kbounds(:,2)) ;
    sWmin = sWmin{1}.*cellConn ;
    sWmin = full(min(sWmin(sWmin>0))) ;
    if isempty(sWmin); sWmin=1 ; end % robustness (case 1 cell)
else
    sWmin = 1;
end

% Various values
Ctr = getTraceConstantL2(assembler) ;
cellSize = getCellSize(assembler) ;
shapeRatio = min(cellSize)^2/prod(cellSize) ;

% Compute penalty lower bound and apply factor
penalty_inf = 2*kPlusMax*Ctr^2*aWmax^2*shapeRatio/(kMinusMin*sWmin) ;
penalty = penalty_inf*factor ;
end

% $$$ function cellConn = mesoConnectivity(cellNum,periodicity)
% $$$ %  cellConn = mesoConnectivity(cellNum,periodicity)
% $$$ cellConn = zeros(prod(cellNum));
% $$$ for direction = 1:2
% $$$     N = cellNum(direction) ;
% $$$     cellConnDir = sparse(1:N-1,2:N,ones(N-1,1),N,N) ;
% $$$     if periodicity(direction)
% $$$         % If periodic boundary condition, connect final cell to first cell
% $$$         cellConnDir(N,1) = 1 ;
% $$$     end
% $$$     orthDir = 1+mod(direction+2,2) ; % 2 if direction=1, 1 if direction=2
% $$$     factors = {eye(cellNum(orthDir)),cellConnDir} ;
% $$$     cellConn = cellConn | kron(factors{direction},factors{orthDir}) ;
% $$$ end
% $$$ end
