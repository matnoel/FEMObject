tsp = 5 ; % time span
tst = 1 ; % number of steps
ns = ceil(tsp/tst) ; % time step
bdf = 1 ; % 1 or 2
useW = true ;
penalty = 1e2 ;

influx = cell(ns,1) ;
for t = 1:ns
    influx{t} = @(x) 1*(x(:,1)<=t/ns) ; % step function
end

model = QPModel('order',2,...
    'cellNum',[20 4],...
    'cellSize',[1 1],...
    'elementSize',[.1 1/3],...
    'tolSVD',1e-4,...
    'verbose',true) ;

diam = (.1:.1:.5)' ;
minK = 1 ;
valK = diam*1e3-minK ;
rK = numel(diam) ;
patterns = struct('name',repmat({'bar'},rK,1),...
    'value',num2cell(valK),...
    'size',mat2cell([diam ones(rK,1)],ones(rK,1),2) ,...
    'center',repmat({[]},rK,1),...
    'offset',repmat({[]},rK,1)) ;
Kmicro = drawCellPattern(getCellCoord(model),patterns)+1 ;
% dist = dealMultinomial([.5 .25 .1 .1 .05],getCellNb(model)) ;
K = distributeMicro(model,dist,Kmicro) ;
% figure; plot(model,K); colorbar
bndK = mesoBounds(model,K) ;
if isempty(penalty)
    penalty = 1.5*penaltyBound(model,bndK,true,true,false) ;
end

ptTol = 1e-6 ;
bc = [QPBC(2,influx{1},@(x) round(x(:,2)/ptTol)*ptTol == 1) ; ...
    QPBC(1,0,@(x) ismember(round(x(:,1)/ptTol)*ptTol,[0 1])) ; ...
    QPBC(2,0,@(x) round(x(:,2)/ptTol)*ptTol == 0) ] ;

pb = QPProblem(model,K,0,bc,'tolerance',1e-1,'useWeights',useW,...
    'penalty',penalty,'verboseUpdater',false) ;

%% Resolutions

sol = cell(1,ns) ;
out = QPProblem.qpOutputStructure(ns) ;
eyeT = TuckerLikeTensor.eye(tensorSize(model)) ;

% Initial solution set to zero
% Loop on time steps
for n = 1:ns
    bc(1).value = QPBC.formatValue(influx{n}) ;
    [swipL,swipR] = swipOperator(model,K,0,bc,penalty,useW,useW) ;
    if bdf==1
        lhs = eyeT+tst*swipL ;
        rhs = tst*swipR ;
        if n>1
            rhs = rhs + sol{n-1} ;
        end
    else % bdf == 2
        if n==1
            lhs = eyeT + .5*tst*swipL ;
            rhs = .5*tst*swipR ;
        else
            lhs = 3*eyeT + 2*tst*swipL ;
            rhs = 2*tst*swipR + 4*sol{n-1} ;
            if n>2
                rhs = rhs - sol{n-2} ;
            end
        end
    end
    if n>1
        pb.initialPoint = sol{n-1} ;
    end
    [sol{n},out(n)] = solve(pb,lhs,rhs) ;
end

%% Post-processing

for n = 1:ns
    figure(n)
    plot(model,sol{n});
    title(sprintf('t_%i',n));
%     caxis manual ;
%     caxis([0 12]);
end