function [Khomo,out] = homogenise(pb,K,stagTol,stdTol,method)
% [Khomo,output] = homogenise(pb,K,stagnationTol,standardDeviationTol,method)
% method : (default is 2)
% _ 1 for low-rank approximation without recycling
% _ 2 for low-rank approximation with recycling
% _ 3 for direct FEM resolution
if nargin < 5
    method = 2 ;
    if nargin < 4
        stdTol = 1e-2 ;
    end
end
verbose = true ;
nTrustSTD = 50 ; % trust STD estimates after this many iterations

% Allocating
iterMax = numel(K) ;
stagnation = zeros(iterMax-1,1) ;
std = zeros(iterMax-1,1) ;
times = zeros(iterMax-1,1) ;
corrK = cell(iterMax,1) ; % corrected conductivity
estK = cell(iterMax,1) ; % estimated conductivity (so far)
if method == 3
    corrector = cell(iterMax,1) ;
    solvOut = QPProblem.feOutputStructure([iterMax 1]) ;
else
    corrector = cell(iterMax,2) ;
    solvOut = QPProblem.qpOutputStructure([iterMax 2]) ;
    coord = getCoord(pb.model) ;
    if method == 2
        microCoord=getCellCoord(pb.model); % for modes rotation
    end
end

% Iterations
flag = 0 ;
for n = 1:iterMax
    ifprint(verbose,'Homogenise %i - ',n)
    pb.K = K{n} ;
    if method == 3
        pb.source = 'correctors' ;
        [corrector{n},solvOut(n)] = pb.solveFE(false,2) ;
        times(n) = sum(solvOut(n).time) ;
        corrK{n} = applyCorrector(solvOut(n).operators.diffOp,...
            corrector{n},getcoord(solvOut(n).model.node)) ;
    else
        % Corrector 1
        pb.source = 'corrector1' ;
        if n>1 && method==2
            pb.initialPoint = corrector{n-1,1} ;
            ifprint(verbose,'Compression from rank %i',pb.initialPoint.space.dim(1))
            tr = Truncator('tolerance',getTolSVD(pb.model),...
                'maxRank',pb.initialPoint.space.dim(1)) ;
            pb.initialPoint = tr.truncate(pb.initialPoint) ;
            ifprint(verbose,' to %i\n',pb.initialPoint.space.dim(1))
        end
        [corrector{n,1},solvOut(n,1),greedy] = pb.solve() ;
        times(n) = solvOut(n,1).time+solvOut(n,1).initTime ;
        % Corrector 2
        pb.source = 'corrector2' ;
        if method==2
            if n>1
                pb.initialPoint = corrector{n-1,2} ;
            else
                try
                    pb.initialPoint = rotateModes(corrector{n,1},microCoord) ;
                catch
                    warning('rotateModes failed')
                end
            end
            ifprint(verbose,'Compression from rank %i',pb.initialPoint.space.dim(1))
            tr = Truncator('tolerance',getTolSVD(pb.model),...
                'maxRank',pb.initialPoint.space.dim(1)) ;
            pb.initialPoint = tr.truncate(pb.initialPoint) ;
            ifprint(verbose,' to %i\n',pb.initialPoint.space.dim(1))
        end
        [corrector{n,2},solvOut(n,2)] = pb.solve(greedy.A) ;
        times(n) = times(n) + ...
            solvOut(n,2).time+solvOut(n,2).initTime ;
        cells = repmat(formatIndex(getOrder(pb.model),...
            getCellNum(pb.model),(1:getCellNb(pb.model))'),1,2) ;
        diffOp = bilinFormOperator(pb.model,[1 1 0],cells,K{n}) ;
        corrK{n}=applyCorrector(diffOp,corrector(n,:),coord);
    end
    ifprint(verbose,'%s',formatDuration(times(n)))
    if n==1
        estK{1} = corrK{1} ;
        ifprint(verbose,'\n')
    else
        estK{n} = (estK{n-1}*(n-1) + corrK{n})/n ;
        stagnation(n-1) = norm(estK{n}-estK{n-1})/norm(estK{n-1}) ;
        ifprint(verbose,' - stagnation %.3g',stagnation(n-1))
        corrKvec = cat(3,corrK{:}) ;
        std(n-1) = sqrt(max(diag(var(corrKvec,0,3))))/sqrt(n) ;
        ifprint(verbose,' - std %.3g\n',std(n-1))
        if n>=nTrustSTD && std(n-1)<stdTol
            disp('Standard deviation criterion satisfied')
            flag = 2 ;
            break
        elseif stagnation(n-1)<stagTol
            disp('Stagnation criterion satisfied')
            flag = 1 ;
            break
        end
    end
end
out = struct('flag',flag,...
    'stagnation',stagnation(1:n-1),...
    'standardDeviation',std(1:n-1),...
    'time',times(1:n),...
    'effectiveConductivity',{estK(1:n)},...
    'correctedConductivity',{corrK(1:n)},...
    'corrector',{corrector(1:n,:)},...
    'conductivity',{K(1:n)},...
    'solverOutput',solvOut(1:n,:));
Khomo = estK{n} ;
end