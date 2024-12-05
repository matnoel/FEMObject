function [solution,output,pb] = solveDD(pb)
% [solution,output,pb] = solveDD(pb)
% Alternate global-local algorithm for resolution in case of domain
% decomposition.

algorithmClock = tic ;

% Solver parameters
compressRHS = true ;
relax = .6 ; % relaxation parameter
errTol = getTolerance(pb) ;
stagTol = getTolerance(pb) ;
maxIter = 10 ;

% Get necessary objects
model = getModel(pb) ;
swipLHS = getLHSOperator(pb) ;
swipRHS = getRHSOperator(pb) ;
patches = getPatches(pb) ;
swipLHSPatches = restrictTensor(patches,swipLHS) ;
patchesLHS = getLHSOperator(patches) ;
patchesRHSOp = getRHSOperator(patches) ;
source = getSource(pb) ;
microOp = getMicroOperators(getOperatorAssembler(pb)) ;
penalty = getPenalty(getOperatorAssembler(pb)) ;

% Compression
compressor = Truncator('tolerance',getTolSVD(pb)) ;
lhs = truncate(compressor,swipLHS + patchesLHS) ;
swipRHS = truncate(compressor,swipRHS) ;
if max(patchesRHSOp{1}.space.dim) > 1 % In case of SWIP patches only
    patchesRHSOp{1} = truncate(compressor,patchesRHSOp{1}) ;
end
patchesRHSOp{2} = truncate(compressor,patchesRHSOp{2}) ;
swipLHSPatches = truncate(compressor,swipLHSPatches) ;

% Initialization
output = QPDiffusionProblem.ddOutputStructure(maxIter) ;
llrhsOpll = norm(swipRHS) ;
errorFun = @(x) norm(swipRHS-swipLHS*x)/llrhsOpll ;
stagnationFun = @(old,new) norm(new-old)/norm(new) ;
U = TuckerLikeTensor.zeros(tensorSize(model));
value = U ;
flux = U ;
output.initTime = toc(algorithmClock) ;

% Iterative alternated resolutions
for iter = 1:maxIter
    iterClock = tic ;    
    
    % Update global
    U_prev = U ;
    rhs = swipRHS + swipLHSPatches*U_prev + patchesRHSOp{1}*value...
        + patchesRHSOp{2}*(-flux) ; % flux is *outward* patches
    if compressRHS
        rhs = truncate(compressor,rhs) ;
    end
    pb = updateGreedySolver(pb,lhs,rhs) ;
%     if iter > 1
%         pb = setInitialPoint(pb,U_prev) ; %TODO: why bug ?
%     end
    
    % Solve global problem
    greedy = getGreedySolver(pb) ;
    [U,gOutput] = solve(greedy) ;
    U = relax*U + (1-relax)*U_prev ;
    
    % Solve local problem
    [w,value,flux,lOutput] = solve(patches,microOp,source,U,errTol,penalty) ;
    
    %/DEBUG
    figure
    subplot(2,1,1)
    plot(model,U);colorbar
    title(sprintf('U_%i',iter))
    subplot(2,1,2)
    plot(getPatchModel(patches,1),subTensor(w,getCellList(patches),':'))
    colorbar
    title(sprintf('w_%i',iter))
    %DEBUG\
    
    % Error indicators
    errGlobal = errorFun(U) ;
    stagGlobal = stagnationFun(U_prev,U) ;
    
    % Store iteration data
    output.error(iter) = errGlobal ;
    output.stagnation(iter) = stagGlobal ;
    output.globalSolution{iter} = U ;
    output.globalOutput{iter} = gOutput ;
    output.localSolution{iter} = w ;
    output.localOutput{iter} = lOutput ;
    output.iterTime(iter) = toc(iterClock) ;
    
    % Check convergence
    residualOk = errGlobal < errTol ;
    stagnationOk = stagGlobal < stagTol ;
    % DEBUG
    fprintf('DD: Iteration %3d - Stagnation %.2d - Error %.2d\n',...
        iter,stagGlobal,errGlobal);
    if residualOk && stagnationOk
        output.flag = 3 ;
        break
    elseif residualOk
        output.flag = 1 ;
%         break
    elseif stagnationOk
        output.flag = 2 ;
        break
    end
end

% Trim output
if iter < maxIter
    output.error = output.error(1:iter) ;
    output.stagnation = output.stagnation(1:iter) ;
    output.globalSolution = output.globalSolution(1:iter) ;
    output.globalOutput = output.globalOutput(1:iter) ;
    output.localSolution = output.localSolution(1:iter) ;
    output.localOutput = output.localOutput(1:iter) ;
    output.iterTime = output.iterTime(1:iter) ;
end

% Assemble final solution
solution = U - restrictTensor(patches,U) + w ;

output.time = toc(algorithmClock) ;
end