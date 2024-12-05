model = QPModel('order',2,...
    'cellNum',[15 14],...
    'cellSize',[1 1],...
    'elementSize',[.1 .1],...
    'tolSVD',1e-6,...
    'verbose',false) ;

patterns = struct('name',{'uniform','rectangle'},...
    'value',{1 99},...
    'size',{[] [.25 .25]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 0 1] ;

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[0.9 0.1]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'constantNullification','full') ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

pb = assembleOperators(pb) ;

%%

% Get greedy and clear memory
greedy = getGreedySolver(pb) ;
tr = Truncator('tolerance',getTolSVD(pb)) ;
order = getOrder(pb) ;
clear pb diffAss KAss model

% Prepare for residual minimisation
A = greedy.A ;
b = greedy.b ;
Amr = tr.truncate(A'*A) ;
bmr = tr.truncate(A'*b) ;
A = tr.truncate(A) ;
b = tr.truncate(b) ;

% Prepare parameters
% #1 minimizeResidual [false,true]
% #2 useDMRG [false,true]
% #3 updateCore [false,true]
% #4 updateTSpace [false,true]
% #5 updateTSpaceDimensions [order,1:order-1,order]
% #6 updater maxIterations [1,2,3,4,5,6,7,8,9,10]
paramDims = [2 2 2 2 3 10] ;
ii = ind2NDsub(paramDims,(1:prod(paramDims))') ;
% Determine some parameters
removed = ii(:,1)==2 | ii(:,4)==1 | ii(:,5)~=3 | ~ismember(ii(:,6),1) ;
ii(removed,:) = [] ;
n = size(ii,1) ;

% Initialisation
iter = 50 ;
rawErrors = zeros(n,iter) ;
rawRanks = zeros(n,iter) ;
rawTimes = zeros(n,iter) ;
for j = 1:iter
fprintf('Set %i - Launching %i simulations\n',j,n)
for i = 1:n
    fprintf('Simulation %i - ',i)
%     greedy.minimizeResidual = ii(i,1)==2 ;
    if ii(i,1)==2 % minimizeResidual
        greedy.A = Amr ;
        greedy.b = bmr ;
    else
        greedy.A = A ;
        greedy.b = b ;
    end
    greedy.solutionUpdater.useDMRG = ii(i,2)==2 ;
    greedy.solutionUpdater.updateCore = ii(i,3)==2 ;
    greedy.solutionUpdater.updateTSpace = ii(i,4)==2 ;
    if ii(i,5)==1
        greedy.solutionUpdater.updateTSpaceDimensions = order ;
    elseif ii(i,5)==2
        greedy.solutionUpdater.updateTSpaceDimensions = 1:order-1 ;
    else
        greedy.solutionUpdater.updateTSpaceDimensions = 1:order ;
    end
    greedy.solutionUpdater.maxIterations = ii(i,6) ;
    [~,out] = solve(greedy) ;
    rawErrors(i,j) = out.error ;
    rawRanks(i,j) = out.iter ;
    rawTimes(i,j) = out.time ;
    fprintf('error %.3g - %i iterations - %s\n',rawErrors(i),rawRanks(i),...
        formatDuration(rawTimes(i)))
end
fprintf('%i simulations finished in %s\n',n,formatDuration(sum(rawTimes(:,j))))
end
fprintf('%i sets of %i simulations finished in %s\n',iter,n,...
    formatDuration(sum(sum(rawTimes)))) 

%%
times = mean(rawTimes,2) ;
ranks = mean(rawRanks,2) ;
errors = max(rawErrors,[],2) ;

[times,ranking] = sort(times) ;
errors = errors(ranking) ;
ranks = ranks(ranking) ;
success = errors < greedy.tolerance & ranks < greedy.maxIterations ;
disp([ranking(success) times(success) ranks(success)])

% Best-performing parameters values : 1 1 1 2 3 1