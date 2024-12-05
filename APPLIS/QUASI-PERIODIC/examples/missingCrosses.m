% Parameters

cellNums = repmat([5 10 15]',1,2) ;
dGfe = false ; % if true, solve dG directly ; if false, use solveFEM.

%% Pre-processing

nSim = size(cellNums,1) ; % number of simulations
times = zeros(nSim,2) ;
timesPP = times ; % pre-processing times (assembling, etc.)
ranks = zeros(nSim,1) ;
errors = ranks ;
sol = cell(nSim,2) ;
output = cell(nSim,2) ;

model = QPModel('order',2,...
    'cellNum',cellNums(1,:),...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','cross'},...
    'value',{1 99},...
    'size',{[] [.2 .2]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 1 0] ;

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[0.9 0.1]) ;
KAss = assemble(KAss) ; % to get conductivity and its bounds

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'useCompression',true,...
    'useAverageWeights',true,...
    'useStabilisationWeights',true,...
    'penalty',[]) ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

%% Processing

for n = 1:nSim
    if n>1
        model = setCellNum(model,cellNums(n,:)) ;
        pb = updateModel(pb,model) ;
    end
    [sol{n,1},output{n,1},pb] = solve(pb,1:getOrder(pb)) ;
    timesPP(n,1) = output{n,1}.initTime ;
    times(n,1) = output{n,1}.time + timesPP(n,1) ;
    ranks(n,1) = output{n,1}.iter ;
    if dGfe
        tic;
        lhsFE = doubleQP(getLHSOperator(pb)) ;
        rhsFE = doubleQP(getRHSOperator(pb)) ;
        timesPP(n,2) = toc ;
        sol{n,2} = lhsFE\rhsFE ;
        times(n,2) = toc;
        errors(n,1) = norm(doubleQP(sol{n,1})-sol{n,2})/norm(sol{n,2}) ;
    else
        [sol{n,2},output{n,2}] = solveFEM(pb,false,2) ;
        timesPP(n,2) = sum(output{n,2}.time(1:2)) ;
        times(n,2) = sum(output{n,2}.time(1:3)) ;
        errors(n,1) = compare2FE(pb,sol{n,1},sol{n,2},...
            @(ref,rel)norm(rel-ref)/norm(ref),'feModel',output{n,2}.model) ;
    end
end

%% Post-processing

femName='FEM';
if dGfe ; femName=['dG-',femName]; end

for n=1:nSim
    fprintf('%i cells: LR %s (rank %i) / %s %s - error %.3g\n',...
        prod(cellNums(n,:)), formatDuration(times(n,1)),ranks(n,1),...
        femName, formatDuration(times(n,2)), errors(n)) ;
end


%% Illustrations

% model2 = QPModel('order',2,...
%     'cellNum',[3 3],...
%     'cellSize',[1 1],...
%     'elementSize',[.01 .01],...
%     'tolSVD',1e-6,...
%     'verbose',true) ;
% patterns = struct('name',{'uniform','cross'},...
%     'value',{1 99},...
%     'size',{[] [.2 .2]} ,...
%     'center',{[] [.5 .5]},...
%     'offset',{[] []}) ;
% patternsTable = [1 1 ; 1 0] ;
% KAss2 = QPConductivityAssembler('model',model2,...
%     'patterns',patterns,...
%     'patternsTable',patternsTable,...
%     'distribution',{setdiff(1:9,5)' 5}) ;
% KAss2 = assemble(KAss2) ;
% diffAss2 = QPDiffusionAssembler('conductivityAssembler',KAss2,...
%     'BC','PBC',...
%     'source','corrector1',...
%     'useCompression',true) ;
% pb2 = QPDiffusionProblem('operatorAssembler',diffAss2,...
%     'tolerance',1e-3);
% 
% % source term
% coord = getDomainCoord(model2) ;
% cSz = getCellSize(model2) ;
% cN = getCellNum(model2) ;
% eSz = getElementSize(model2) ;
% X = reshape(coord(:,1),cN(1)*(1+cSz(1)/eSz(1)),[]) ;
% Y = reshape(coord(:,2),[],cN(2)*(1+cSz(2)/eSz(2))) ;
% assert(all(size(X))==all(size(Y)),'Inconsistent sizes')
% K = doubleQP(getConductivity(KAss2)) ;
% K1 = reshape(K,size(X,1),size(X,2)) ;
% K2 = zeros(size(K1)) ;
% src = divergence(X,Y,K1,K2) ;
% figure
% plot(model2,src(:))
% 
% % approximation
% [sol2,~,pb2] = solve(pb2,1:getOrder(model2)) ;
% figure
% plot(model2,sol2)
% 
% % solution FEM
% solFEM = solveFEM(pb2);
% figure
% plot(model2,sol2)