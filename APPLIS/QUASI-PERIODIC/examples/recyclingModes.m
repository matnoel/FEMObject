%% Pre-processing

refModel = QPModel('order',2,...
    'cellNum',[50 20],...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','rectangle'},...
    'value',{1 99},...
    'size',{[] [.25 .25]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
refKAss = QPConductivityAssembler('model',refModel,...
    'patterns',patterns,...
    'patternsTable',[1 1 ; 1 0],...
    'probability',[0.5 0.5]) ;
refDiffAss = QPDiffusionAssembler('conductivityAssembler',refKAss,...
    'BC','PBC',...
    'source','corrector1',...
    'constantNullification','full') ;

refPb = QPDiffusionProblem('operatorAssembler',refDiffAss,...
    'tolerance',1e-3);

%% Resolution

N = 50 ;
[~,outRM,refPb] = multiSolve(refPb,N) ;
%[~,outFE] = solveFEM(refPb) ;

%% Post-processing

times = zeros(N,2) ;
ranks = zeros(N,1) ;
for i = 1:N
    times(i,1) = outRM(i).time + outRM(i).initTime ;
    ranks(i) = outRM(i).iter ;
end
%times(:,2) = repmat(sum(outFE.time),N,1) ;

figure
hold on
plot(1:N,times(:,1),'-xb')
%plot(1:N,times(:,2),'-or')
ylabel('Computation time (s)')
%legend({'LR','FEM'})

figure
hold on
plot(1:N,cumsum(times(:,1)),'-xb')
%plot(1:N,cumsum(times(:,2)),'-or')
ylabel('Cumulated computation time (s)')
%legend({'LR','FEM'})