N = 30 ;
tol = 10.^linspace(-1,-10,N);

model = QPModel('order',2,...
    'cellNum',[20 20],...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',min(1e-6,tol(end)/100),...
    'verbose',true) ;

patterns = struct('name',{'uniform','rectangle'},...
    'value',{1 99},...
    'size',{[] [.25 .25]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 1 0] ;
Kmicro = drawCellPattern(getCellCoord(model),patterns)*patternsTable ;
Kdist = dealMultinomial([0.9 0.1],getCellNb(model)) ;
K = distributeMicro(model,Kdist,Kmicro) ;

pb = QPProblem(model,K,'corrector1',QPBC(4),'tolerance',tol(1),...
    'penalty',100,'verboseUpdater',false) ;

%% Resolution

ranks = zeros(N,1) ;
errors = zeros(N,1) ;

% Initial run
fprintf('#%i/%i - Tolerance %.3g\n',1,N,tol(1))
[sol,out,greedy] = pb.solve() ;
ranks(1) = sol.space.dim(1) ;
errors(1) = out.error ;

% Iterations
for n=2:N
    fprintf('#%i/%i - Tolerance %.3g\n',n,N,tol(n))
    pb.initialPoint = sol ;
    pb.tolerance = tol(n) ;
    [sol,out,greedy] = pb.solve(greedy.A,greedy.b) ;
    ranks(n) = sol.space.dim(1) ;
    errors(n) = out.error ;
end

%% Post-processing

% matlab2tikz additional arguments
m2tikzArgin = {'showInfo',false,'showwarnings',false,'noSize',true,...
    'extraAxisOptions',{'ymajorgrids=true','grid style ={dotted}'}} ;

figure
plot(ranks,log10(errors),'-db')
xlabel('Approximation''s rank')
ylabel('Precision')
cleanfigure('pruneText',false);
tikzName = 'precisionWRTRank.tikz' ;
% matlab2tikz(tikzName,m2tikzArgin{:})