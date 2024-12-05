maxIter = 7 ;
pow = [(2*(maxIter-1)):-1:(maxIter-1) ; 0:(maxIter-1)]' ;

model2 = QPModel('order',2,...
    'cellNum',[2 2].^pow(1,:),...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','rectangle'},...
    'value',{1 99},...
    'size',{[] [.25 .25]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 1 0] ;

KAss = QPConductivityAssembler('model',model2,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[0.9 0.1]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'constantNullification','full',...
    'useAverageWeights',false,...
    'useStabilisationWeights',false) ;

pb2 = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-2);

times = zeros(maxIter,2) ;
ranks = zeros(maxIter,2) ;
for i = 1:maxIter
    % Order 2
    tic ;
    model2 = setCellNum(model2,[2 2].^pow(i,:)) ;
    pb2 = updateModel(pb2,model2) ;
    updateTime = toc ;
    [sol2,out2] = solve(pb2) ;
    times(i,1) = out2.time + out2.initTime + updateTime;
    ranks(i,1) = out2.iter ;
    % Order 3
    tic ;
    model3 = setOrder(model2,3) ;
    pb3 = updateModel(pb2,model3) ;
    updateTime = toc ;
    [sol3,out3] = solve(pb3) ;
    times(i,2) = out3.time + out3.initTime + updateTime;
    ranks(i,2) = out3.iter ;
end

%%

ratio = 2*prod(2.^pow,2)./sum(2.^pow,2) ;
figure ; hold on
plot(ratio,times(:,1),'-+b')
plot(ratio,times(:,2),'-+r')
legend('order2','order3')
xlabel('N_1 + N_2 / N_1 \times N_2')
ylabel('Time (s)')

figure ; hold on
plot(ratio,ranks(:,1),'-+b')
plot(ratio,ranks(:,2),'-+r')
legend('order2','order3')
xlabel('2\times N_1\times N_2/(N_1+N_2)')
ylabel('Rank')
