function [sol,output] = solve(pb,localSolvers,tolerance,reference,relaxation,fullModel)
% [sol,output] = solve(pb,localSolvers,tolerance,reference,relaxation,fullModel)

if nargin < 6
    fullModel = [] ;
    if nargin < 5
        relaxation = [] ;
        if nargin < 4
            reference = [] ;
        end
    end
end

minRelaxation = 1e-2 ; % hard coded lower bound on allowed values
maxRelaxation = 1e1 ; % hard coded upper bound on allowed values

aitken = isempty(relaxation) || ~isnumeric(relaxation) || ...
    ~isscalar(relaxation) || relaxation<=0 || relaxation>=2 ;
patchNb = numel(pb.patch) ;

% Process criterion input
stagnationCriterion = isempty(reference) || numel(tolerance)>1 ;
stagTolerance = tolerance ;
referenceCriterion = ~isempty(reference) ;
if referenceCriterion
    referenceNorm = norm(reference) ;
    refTolerance = tolerance(1) ;
    if numel(tolerance)>1
        stagTolerance = tolerance(2) ;
    end
end

%% Compute various patch's boundary coordinates and transfer operators

% Build *local* boundary coordinates
% and *global* internal boundary coordinates
lBoundaryCoord = cell(patchNb,1) ;
gIntBoundaryCoord = cell(patchNb,1) ;
for p = 1:patchNb
    % Local coordinates
    lBoundaryCoord{p} = getDomainCoord(pb.patch(p).model) ;
    pSz = max(lBoundaryCoord{p},[],1) ; % patch size
    locBoundary = lBoundaryCoord{p}(:,2)<pb.coordTol | ... % bottom edge
        pSz(1)-lBoundaryCoord{p}(:,1)<pb.coordTol | ... % right edge
        pSz(2)-lBoundaryCoord{p}(:,2)<pb.coordTol | ... % top edge
        lBoundaryCoord{p}(:,1)<pb.coordTol ; % left edge
    lBoundaryCoord{p} = lBoundaryCoord{p}(locBoundary,:) ;
    % Global coordinates (=local+offset)
    offset = min(pb.patchCells{p}) ; % index of lewer left patch cell
    offset = formatIndex(3,pb.model.cellNum,offset) ; % as subscripts
    offset = (offset-1)*diag(getCellSize(pb.model)) ; % then scale
    gIntBoundaryCoord{p} = [lBoundaryCoord{p}(:,1)+offset(1) , ...
        lBoundaryCoord{p}(:,2)+offset(2)] ;
end

% Make two list of cells of patches' external boundary
% One by (direction,patch} and one by direction.
patchExtCells = cell(4,patchNb) ;
extCells = cell(4,1) ; % list by direction
cCoord = getCellCoord(pb.model) ;
cSz = max(cCoord,[],1) ;
locCellBoundary = {cSz(1)-cCoord(:,1)<pb.coordTol ... % right edge
    cSz(2)-cCoord(:,2)<pb.coordTol ... % top edge
    cCoord(:,1)<pb.coordTol ... % left edge
    cCoord(:,2)<pb.coordTol} ; % bottom edge
for p = 1:patchNb
    patchExtCells(:,p) = boundaryCells(pb.model,pb.patchCells{p},true) ;
    extCells = cellfun(@union,extCells,patchExtCells(:,p),'UniformOutput',false);
end

% Build *global* external boundary coordinates
coordFull = getDomainCoord(pb.model) ;
gExtBoundaryCoord = repmat({zeros([0 2])},patchNb,1) ;
for p = 1:patchNb
    for d=1:4
        new = mesoIndicatorT(pb.model,patchExtCells{d,p}) ;
        new.space.spaces{end} = double(locCellBoundary{d}) ;
        gExtBoundaryCoord{p} = gExtBoundaryCoord{p} + new ;
    end
    loc = doubleQP(gExtBoundaryCoord{p})==1 ;
    gExtBoundaryCoord{p} = full(coordFull) ;
    gExtBoundaryCoord{p}(~loc,:) = -1 ;
    % Extend to whole domain with zeros and correct order
%     gExtBoundaryCoord{p} = relativeSort(gExtBoundaryCoord{p},...
%         round(gExtBoundaryCoord{p}/pb.coordTol)*pb.coordTol,coordFull,true) ;
end

% Build transfer operators on patch boundaries
% From exterior side in global discontinuous space to interior side in
% local trace space, and conversely
ext2IntOp = cell(patchNb,1) ;
int2ExtOp = cell(patchNb,1) ;
for p = 1:patchNb
    % Transfer operator from exterior to interior
    % Round coordinates to tolerance for comparison
    extCoord = round(gExtBoundaryCoord{p}/pb.coordTol)*pb.coordTol ;
    intCoord = round(gIntBoundaryCoord{p}/pb.coordTol)*pb.coordTol ;
    % Compare. intCoord(matchedInt,:)=extCoord(ext2Int,:)
    [matchedInt,ext2Int] = ismember(intCoord,extCoord,'rows') ;
    if ~all(matchedInt)
        warning('Not all internal boundary nodes of patch %i were matched',p);
    end
    % Matching meshes assumed, so nodes not listed in ext2Int are
    % duplicates. List them and find matching internal nodes.
    duplicateExt = setdiff((1:size(extCoord,1))',union(unique(ext2Int),...
        find(all(extCoord==-1,2)))) ;
    [matchedExt,int2Ext] = ismember(extCoord,intCoord,'rows') ; % needed for int2ExtOp
    if any(any(extCoord(~matchedExt,:)>=0))
        warning('Not all external boundary nodes of patch %i were matched',p);
    end
    % Build matrix such that ext2IntOp{p}(m,n) = 1 if 
    % extCoord(n,:)==intCoord(m,:) and else 0.
    i = [find(matchedInt) ; int2Ext(duplicateExt)] ;
    j = [ext2Int ; duplicateExt] ;
    ext2IntOp{p} = sparse(i,j,1,size(intCoord,1),size(extCoord,1)) ;
    
    % Transfer operator from interior to exterior
%     duplicateInt = setdiff((1:size(intCoord,1))',unique(int2Ext)) ;
%     assert(isempty(duplicateInt),'I found duplicate interior nodes') ; % should not happen
%     int2ExtOp{p} = sparse(find(matchedExt),int2Ext(int2Ext>0),1,...
%         size(extCoord,1),size(intCoord,1)) ;
    int2ExtOp{p} = ext2IntOp{p}' ;
    
    % Normalise so the sum of each line be one (i.e. smooth external nodes'
    % values by averaging) in ext2IntOp only
    ext2IntOp{p} = diag(sum(ext2IntOp{p},2).^-1)*ext2IntOp{p} ;
end

%TODO: Transfer BC on patches from global to local

%% Assemble operators

% Penalisation
% (To ensure the same is used for global and patch operators)
periodicity = [isPeriodic(pb.bc,1) isPeriodic(pb.bc,2)] ;
penalty = pb.penalty ;
if isempty(penalty)
    penalty = 1.5*penaltyBound(pb.model,mesoBounds(pb.model,pb.K),...
        false,false,periodicity) ;
end
for p = 1:patchNb
    for ibc = 1:numel(pb.patch(p).bc)
        if pb.patch(p).bc(ibc).type==1 && isempty(pb.patch(p).bc(ibc).factor)
            pb.patch(p).bc(ibc).factor = 100*penalty ;
        end
    end
end

% SWIP global operators
[globalLHS,swipRHS]=swipOperator(pb.model,pb.K,pb.source,pb.bc,penalty,...
    pb.useWeights,pb.useWeights) ;
% Set swipRHS to zero on patches (affects only source operator and
% constant nullification operator, unless patch intersects boundary)
patchCellsList = pb.allPatchCells ;
nonPatchCells = setdiff(1:getCellNb(pb.model),patchCellsList) ;
swipRHS = restrictTensor(pb.model,swipRHS,nonPatchCells,false) ;

% SWIP patch operator
cNbG = getCellNb(pb.model) ;
connecP = {sparse(cNbG,cNbG) sparse(cNbG,cNbG)} ; % patches connectivity
for d=1:2
    connecd = mesoConnectivity(pb.model,d,periodicity(d)) ;
    connecP{d}(patchCellsList,:) = connecd(patchCellsList,:) ;
    connecP{d}(:,patchCellsList) = connecd(:,patchCellsList) ;
end
consOpPatch = consistencyOperator(pb.model,pb.K,connecP,[]) ;
patchSWIPOp = bilinFormOperator(pb.model,[1 1 0],patchCellsList,pb.K) ...
    + penalty*stabilisationOperator(pb.model,1,connecP,[]) ...
    - consOpPatch - consOpPatch' ;
if undefinedConstant(pb.bc)
    patchSWIPOp = patchSWIPOp + bilinFormOperator(pb.model,0,patchCellsList,1) ;
end

% Inner product operator on patches' exterior boundaries
gBoundaryOp = bilinFormOperator(pb.model,0,extCells,1) ;

%% Global-Local loop

tr = Truncator('tolerance',getTolSVD(pb.model)) ;
globalLHS = tr.truncate(globalLHS) ;
pb.compress = false ;

% Allocate output structure
output.flag = -1 ;
output.time = zeros(pb.maxIterationsGL,1+patchNb) ;
output.relaxation = zeros(pb.maxIterationsGL,1) ;
output.error = zeros(pb.maxIterationsGL,1) ;
output.stagnation = zeros(pb.maxIterationsGL,1) ;
output.globalSolution = cell(pb.maxIterationsGL,1) ;
output.globalOutput = QPProblem.qpOutputStructure([pb.maxIterationsGL 1]) ;
output.localSolution = cell(pb.maxIterationsGL,patchNb) ;
output.localOutput = QPProblem.feOutputStructure([pb.maxIterationsGL patchNb]) ;
% only FE solver implemented

% Initialise then loop
prevGSol = TuckerLikeTensor.zeros(pb.model.tensorSize) ;
globalSol = TuckerLikeTensor.zeros(pb.model.tensorSize) ;
flux = TuckerLikeTensor.zeros(pb.model.tensorSize) ;
localSol = cell(patchNb,1) ;
for n = 1:pb.maxIterationsGL
    % % Global step
    gClock = tic ;
    % Update global problem
    prev2GSol = prevGSol ;
    prevGSol = globalSol ;
    globalRHS = swipRHS + patchSWIPOp*prevGSol - gBoundaryOp*flux ;
    tr.maxRank = max(globalRHS.space.dim) ;
    globalRHS = tr.truncate(globalRHS) ;
    if n>1
        tr.maxRank = max(prevGSol.space.dim) ;
        pb.initialPoint = tr.truncate(prevGSol) ;
    end
    % Solve global problem
    prevCompress = pb.compress ;
    pb.compress = false ;
    [globalSol,output.globalOutput(n)] = solve@QPProblem(pb,globalLHS,globalRHS) ;
    pb.compress = prevCompress ;
            
    % % Relaxation step
    if aitken
        if n>3
            prevDelta = delta ;
            delta = globalSol-prevGSol ;
            deltaDiff = delta-prevDelta ;
            relaxation = -relaxation*dot(deltaDiff,prevDelta)/dot(deltaDiff,deltaDiff) ;
        elseif n==3 % Aitken initialisation
            prevDelta = delta ;
            delta = globalSol-prevGSol ;
            deltaDiff = delta-prevDelta ;
            relaxation = -dot(deltaDiff,prevGSol-prev2GSol)/dot(deltaDiff,deltaDiff) ;
        else % not enough iterations for Aitken
            relaxation = 1;
            delta = globalSol-prevGSol ;
        end
    end
    % Project onto admissible interval
    relaxation = min(maxRelaxation,max(minRelaxation,relaxation)) ;
    % Relax
    globalSol = relaxation*globalSol + (1-relaxation)*prevGSol ;
    output.relaxation(n) = relaxation ;
    output.globalSolution{n} = globalSol ;
    output.time(n,1) = toc(gClock) ;
    
    % % Local step
    for p = 1:patchNb
        lClock = tic ;
        % Update local problem BC (3 steps)
        % 1_First extract value from globalSol on external boundary
        patchBoundaryInd = [] ;
        for d=1:4 % must be done in same order as for gExtBoundary to use ext2IntOp
            new = mesoIndicatorT(pb.model,patchExtCells{d,p}) ;
            new.space.spaces{end} = double(locCellBoundary{d}) ;
            patchBoundaryInd = patchBoundaryInd + new ;
        end
        gSolBoundary = globalSol.*patchBoundaryInd ;
%         gSolBoundary = relativeSort(gSolBoundary,...% Extend to whole domain
%             gExtBoundaryCoord{p},coordFull,pb.coordTol,true) ; % with zeros and correct order
        % 2_Then project it on internal boundary space (with averaging)
        gSolBoundary = ext2IntOp{p}*doubleQP(gSolBoundary) ;
        % 3_Finally update BC value
        relativeCoord = lBoundaryCoord{p}*diag(max(lBoundaryCoord{p},[],1).^-1) ;
        pb.patch(p).bc(1).value = @(x) relativeSort(gSolBoundary,...
            relativeCoord,x,pb.coordTol,true) ;
        
        % Now solve
        if localSolvers(p) == 1
            [localSol{p},output.localOutput(n,p)] = solveFE(pb.patch(p)) ;
        else
            error('Requested local solver not implemented')
        end
        output.time(n,1+p) = toc(lClock) ;
    end
    output.localSolution(n,:) = localSol' ;
    lClock = tic ;
    if n==1
        [flux,lBoundaryOp,lBoundaryFluxOp] = patchFlux(pb,localSol,[],[],int2ExtOp) ;
    else
        flux = patchFlux(pb,localSol,lBoundaryOp,lBoundaryFluxOp,int2ExtOp) ;
    end
    % Account for flux computation time (assume equal time per patch)
    output.time(n,2:end) = output.time(n,2:end) + toc(lClock)/patchNb ;
    
    % % Stopping criterion
    % Stagnation
    prevGSolNorm = norm(prevGSol) ;
    if prevGSolNorm>0
        stagnation = norm(globalSol-prevGSol)/prevGSolNorm ;
    else
        stagnation = NaN ;
    end
    ifprint(pb.verboseGL,'#%i GL - stag %.2g',n,stagnation)
    if stagnationCriterion && stagnation < stagTolerance
        output.flag = 1 ;
    end
    output.stagnation(n) = stagnation ;
    % Reference
    if referenceCriterion
        % Assemble complete solution
        sol = assembleGlobalLocal(pb,globalSol,localSol,fullModel) ;
        refError = norm(sol-reference)/referenceNorm ;
        ifprint(pb.verboseGL,' - error %.2g',refError)
        if refError < refTolerance
            output.flag = 2 ;
        end
        output.error(n) = refError ;
    end
    ifprint(pb.verboseGL,'\n')
    
    % % DEBUG
    if mod(n,4)==-1
%         global refSolG refSolL refSol refLambdaL
        viewAngle = [0 0 1] ;
        figure
        subplot(2,2,1)
        plot(pb.model,globalSol) ; colorbar ; view(viewAngle) ;
        title(sprintf('Global solution at step %i',n))
        for j = 1:3
            subplot(2,2,1+j)
            plot(pb.patch(j).model,localSol{j}) ; colorbar ; view(viewAngle) ;
            title(sprintf('Local solution %i',j))
        end
%         subplot(2,2,3)
%         lFlux = ext2IntOp{p}*doubleQP(flux.*patchBoundaryInd) ;
%         lFlux = relativeSort(lFlux,lBoundaryCoord{p},...
%             getDomainCoord(pb.patch(p).model),pb.coordTol,true) ;
%         plot(pb.patch(p).model,lFlux) ; colorbar ; %view(viewAngle) ;
%         title(sprintf('Local flux %i',n))
%         subplot(2,2,4)
%         plot(pb.model,doubleQP(sol)) ; colorbar ; view(viewAngle) ;
%         title(sprintf('Assembled solution %i',n))
    end
    % %
    
    % Stop
    if output.flag>0
        break
    end
end

%% Process results

% Assemble complete solution
if ~referenceCriterion
    sol = assembleGlobalLocal(pb,globalSol,localSol,fullModel) ;
end

% Trim outputs
output.time = output.time(1:n,:) ;
output.relaxation = output.relaxation(1:n) ;
output.error = output.error(1:n) ;
output.stagnation = output.stagnation(1:n) ;
output.globalSolution = output.globalSolution(1:n) ;
output.globalOutput = output.globalOutput(1:n) ;
output.localSolution = output.localSolution(1:n,:) ;
output.localOutput = output.localOutput(1:n,:) ;
end

function [flux,M,N,P] = patchFlux(pb,w,M,N,P)

if nargin < 5
    P = [] ;
    if nargin < 4
        N = [] ;
        if nargin < 3
            M = [] ;
        end
    end
end

buildM = isempty(M) ;
buildN = isempty(N) ;
buildP = isempty(P) ;

patchNb = numel(pb.patch) ;
flux = [] ;
if buildM
    M = cell(patchNb,1) ;
end
if buildN
    N = cell(patchNb,1) ;
end
if buildP
    P = cell(patchNb,1) ;
end
for n = 1:patchNb
    if true%buildM
        M{n} = bilinFormOperator(pb.patch(n).model,0,boundaryCells(pb.patch(n).model)) ;
        M{n} = doubleQP(M{n}) ;
        % Reduce space to trace
        [locBdry,~]=find(M{n}); % boundary nodes
        locBdry = unique(locBdry) ;
        M{n} = M{n}(locBdry,:) ;
        M{n} = M{n}(:,locBdry) ;
    end
    if buildN
        N{n} = bilinFormOperator(pb.patch(n).model,[0 1 0],...
            boundaryCells(pb.patch(n).model),pb.patch(n).K) ;
        N{n} = doubleQP(N{n}) ;
        % Reduce *destination* space to trace
        if ~buildM % else, rows are the same
            [locBdry,~]=find(N{n});
            locBdry = unique(locBdry) ;
        end
        N{n} = N{n}(locBdry,:) ;
    end
    if buildP
        %TODO
        error('Not implemented. patchFlux expected transfer operators')
    end
    % Compute flux
    w{n} = sunder(pb.patch(n).model,w{n}) ; 
    lFlux = M{n}\(N{n}*w{n}) ;
%     % Transfer from local to global and from internal to external boundary
%     patchLoc = untensorize(pb.model,mesoIndicatorT(pb,pb.patchCells{p}))>0 ;
%     currentFlux = zeros(size(patchLoc)) ;
%     currentFlux(patchLoc) = lFlux ;
%     currentFlux = sunder(pb.model,currentFlux) ;
    currentFlux = P{n}*lFlux ;

    % Tensorize (no compression) then sum
    flux = flux + tensorize(pb.model,currentFlux,0) ;
end

% Lossless compression
flux = orth(flux) ;

end

% Legacy (implemented separately)
% function completeSol = assembleGlobalLocal(pb,globalSol,localSols,fullModel)
% % completeSol = assembleGlobalLocal(pb,globalSol,localSols,fullModel)
% % _pb: QPDDProblem
% % _globalSol: TuckerLikeTensor
% % _localSols: cell array of double arrays
% % _fullModel: MODEL
% % _completeSol: TuckerLikeTensor if fullModel is empty, else double array
% 
% if nargin<4
%     fullModel=[];
% end
% 
% % Set global solution to zero on patches
% nonPatchCells = setdiff(1:getCellNb(pb.model),pb.allPatchCells) ;
% globalRestricted = restrictTensor(pb.model,globalSol,nonPatchCells,...
%     pb.compress) ;
% 
% if isempty(fullModel) % Then global model is complete model. Process as tensors.
%     % Extend local solutions with zeros, as tensors
%     localExtended = [] ;
%     for i = 1:numel(localSols)
%         localSols{i} = tensorize(pb.patch(i).model,localSols{i}) ;
%         localExtended = localExtended + ...
%             extendTensor(pb.model,localSols{i},pb.patchCells{i}) ;
%     end
% else % Then fullModel is complete model (class MODEL). Process as double
%     % Get global solution outside patches
%     patchLoc = doubleQP(mesoIndicatorT(model,pb.allPatchCells)) ;
%     globalRestricted = doubleQP(globalSol) ;
%     globalRestricted = globalRestricted(~patchLoc) ;
%     % Reorder to fullModel node numbering (node assumed to match outside
%     % patches), complete with zeros
%     gCoord = getDomainCoord(pb.model) ; % global coord
%     gCoord = gCoord(~patchLoc) ; % restrict to match globalRestricted
%     fCoord = getcoord(getnode(fullModel)) ; % full coord
%     globalRestricted = relativeSort(globalRestricted,gCoord,fCoord,...
%         pb.coordTol,true) ;
%     % Reorder local solutions to fullModel node numbering (node assumed to
%     % match inside patches), complete with zeros
%     localExtended = zeros(size(globalRestricted)) ;
%     for i = 1:numel(localSols)
%         % Compute local coordinate offset with respect to global domain
%         offset = min(pb.patchCells{i}) ; % index of lewer left patch cell
%         offset = formatIndex(3,pb.model.cellNum,offset) ; % as subscripts
%         offset = (offset-1)*diag(getCellSize(pb.model)) ; % then scale
%         lCoord = offset + getDomainCoord(pb.patch(i).model);
%         % Reorder (assuming match) and complete with zeros
%         localExtended = localExtended + relativeSort(localSols{i},...
%             lCoord,fCoord,pb.coordTol,true) ;
%     end
% end
% 
% % Final assembling
% completeSol = globalRestricted + localExtended ;
% % TuckerLikeTensor if fullModel is empty, else double array
% end