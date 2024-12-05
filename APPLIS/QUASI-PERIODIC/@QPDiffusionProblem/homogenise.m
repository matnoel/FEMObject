function [Khomo,output] = homogenise(pb,solverChoice,iterMax,stagnationTol,K)
% [Khomo,output] = homogenise(pb,solverChoice,iterMax,stagnationTol,K)

if nargin < 5
    K = [] ;
    if nargin < 4
        stagnationTol = getTolerance(pb) ;
        if nargin < 3
            iterMax = sqrt(getCellNb(pb)) ;
        end
    end
end

if isempty(K)
    K = cell(iterMax,1) ; 
    for i = 2:iterMax % not used on first iteration
        K{i} = getConductivity(updateDistribution(pb)) ;
    end
end

model = getModel(pb) ;
tr = Truncator('tolerance',getTolSVD(model)) ;

% initial run and pre-allocating
stagnations = zeros(iterMax-1,1) ;
corrector = cell(iterMax,2) ;
correctedK = cell(iterMax,1) ;
effectiveK = cell(iterMax,1) ;
fprintf('Homogenise 1 - ') ;
if solverChoice==1
%     solverOutput = QPDiffusionProblem.qpOutputStructure([iterMax 2]) ;
    solverOutput = GreedyLinearSolver.outputStructure([iterMax 2]) ;
    corrPb = setSource(pb,'corrector1') ;
    corrPb = updateConductivity(corrPb,K{1}) ;
    [lhs,rhs] = swipOperator(model,K{1},'corrector1',true,[],false,false) ;
    lhs = lhs+bilinFormOperator(model,0) ;
    corrPb = updateGreedySolver(corrPb,lhs,rhs) ;
    [corrector{1,1},solverOutput(1,1)] = solve(getGreedySolver(corrPb)) ;
%     [corrector{1,1},solverOutput(1,1),corrPb] = solve(corrPb) ;
%    corrPb = updateSource(pb,'corrector2') ;
%    [corrector{1,2},solverOutput(1,2)] = solve(corrPb) ;
    corrector{1,2} = TuckerLikeTensor.zeros(tensorSize(getModel(pb))) ;
    correctedK{1} = applyCorrector(corrPb,corrector(1,:)) ;
else
    solverOutput = QPDiffusionProblem.feOutputStructure([iterMax 1]) ;
    homogPb = setSource(pb,'correctors') ;
    [corr,solverOutput(1,1)] = solveFEM(homogPb,false,2) ;
    corrector(1,:) = mat2cell(corr,size(corr,1),ones(1,size(corr,2))) ;
    correctedK{1} = applyCorrectorFE(getcoord(getnode(solverOutput(1,1).model)), ...
        solverOutput(1,1).operators.diffOp,corrector(1,:)) ;
end
effectiveK{1} = correctedK{1} ;
fprintf('stagnation NaN\n') ;

% Iterative runs
flag = 0 ;
for iter = 2:iterMax

    fprintf('Homogenise %i - ',iter)
    if solverChoice==1
        % Corrector problem 1
        corrPb = setSource(pb,'corrector1') ;
        corrPb = updateConductivity(corrPb,K{iter}) ;
        [lhs,rhs] = swipOperator(model,K{iter},'corrector1',true,[],false,false) ;
        lhs = lhs+bilinFormOperator(model,0) ;
        corrPb = updateGreedySolver(corrPb,tr.truncate(lhs),tr.truncate(rhs)) ;
        [corrector{iter,1},solverOutput(iter,1)] = solve(getGreedySolver(corrPb)) ;
        %        corrPb = setInitialPoint(corrPb,corrector{iter-1,1}) ;
%         [corrector{iter,1},solverOutput(iter,1),corrPb] = solve(corrPb) ;
        
        % Corrector problem 2
%         [corrPb,assemblerTime] = updateSource(corrPb,'corrector2') ;
%         %        corrPb = setInitialPoint(corrPb,corrector{iter-1,2}) ;
%         [corrector{iter,2},solverOutput(iter,2)] = solve(corrPb,[]) ;
%         solverOutput(iter,2).initTime = assemblerTime + solverOutput(iter,2).initTime ;
        corrector{iter,2} = TuckerLikeTensor.zeros(tensorSize(getModel(pb))) ;
        correctedK{iter} = applyCorrector(corrPb,corrector(iter,:)) ;
    
    else % Solve both at once
        homogPb = setConductivity(homogPb,K{iter}) ; 
        [corr,solverOutput(iter,1)] = solveFEM(homogPb,false,2) ;
        corrector(iter,:) = mat2cell(corr,size(corr,1),ones(1,size(corr,2))) ;
        correctedK{iter} = applyCorrectorFE(getcoord(getnode(solverOutput(iter,1).model)), ...
            solverOutput(iter,1).operators.diffOp,corrector(iter,:)) ;
    end
    
    % Correction, storage and stagnation
    effectiveK{iter} = (effectiveK{iter-1}*(iter-1) + correctedK{iter})/iter ;
    stagnations(iter-1) = norm(effectiveK{iter}-effectiveK{iter-1})/norm(effectiveK{iter-1}) ;
    fprintf('stagnation %.3g\n',stagnations(iter-1))
    if stagnations(iter-1) < stagnationTol
        flag = 1 ;
        break
    end
end

if solverChoice==1
%     times = sum([solverOutput(1:iter,1).initTime ; solverOutput(1:iter,1).time ; ...
%         solverOutput(1:iter,2).initTime ; solverOutput(1:iter,2).time])' ;
    times = sum([cat(1,solverOutput(1:iter,1).time) ...
        cat(1,solverOutput(1:iter,2).time)]) ;
else
    times = cat(1,solverOutput(1:iter,1).time) ;
    times = sum(times,2) ;
end

output = struct('flag',flag,...
    'stagnations',stagnations(1:iter-1),...
    'times',times,...
    'effectiveConductivities',{effectiveK(1:iter)},...
    'correctedConductivites',{correctedK(1:iter)},...
    'correctors',{corrector(1:iter,:)},...
    'conductivities',{K(1:iter)},...
    'solverOutputs',solverOutput(1:iter,:));
Khomo = effectiveK{end} ;
end

function corrected = applyCorrectorFE(coord,operator,corrector)
% Apply corrector fields to get (apparent) homogenised quantity

corrected = zeros(2) ;
for i=1:2
    for j=1:2
        corrected(i,j) = (coord(:,i)+corrector{i})'*operator*coord(:,j) ;
    end
end

% scaling
domainMeasure = prod(max(coord)-min(coord)) ;
corrected = corrected/domainMeasure ;
% TODO : why is scaling needed ?
end