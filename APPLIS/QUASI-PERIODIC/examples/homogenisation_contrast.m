% Parameters
caseChoice = 1 ; % 1 unidirectional | 2 bidirectional
maxIter = 20 ;
stagnationTol = 0 ;
approxTol = 1e-2 ;
mesoSize = 400 ;
contrastK = [2 1e2 1e4] ;
proba = .1 ;

%% Pre-processing

switch caseChoice
    case 1
        model = QPModel('order',2,...
            'cellNum',[mesoSize(1) 1],...
            'cellSize',[1 mesoSize(1)],...
            'elementSize',[1 mesoSize(1)]/20,...
            'tolSVD',1e-6,...
            'verbose',false) ;
        patterns = struct('name',{'uniform','bar'},...
            'value',{1 contrastK(1)-1},...
            'size',{[] [.5 1]} ,...
            'center',{[] [.5 .5]},...
            'offset',{[] []}) ;
        patternsTable = [1 1 ; 1 0] ;
        
    case 2
        model = QPModel('order',2,...
            'cellNum',floor(sqrt(mesoSize(1)*[1 1])),...
            'cellSize',[1 1],...
            'elementSize',[.05 .05],...
            'tolSVD',1e-6,...
            'verbose',false) ;
        patterns = struct('name',{'uniform','rectangle'},...
            'value',{1 contrastK(1)-1},...
            'size',{[] sqrt([.5 .5])} ,...
            'center',{[] [.5 .5]},...
            'offset',{[] []}) ;
        patternsTable = [1 1 ; 1 0] ;
end

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[1-proba(1) proba(1)]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1') ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,'tolerance',approxTol(1)) ;

nParam = numel(contrastK) ;
conductivities = cell(maxIter,nParam) ;
outQP = cell(1,nParam) ;
outFE = cell(1,nParam) ;
dist = cell(1,maxIter) ;
cellCoord = getCellCoord(model) ;
for j = 1:nParam
    patterns(2).value = (contrastK(j)-1)*patterns(1).value ;
    microPatterns = drawCellPattern(cellCoord,patterns)*patternsTable ;
    for i = 1:maxIter
        if j==1
            dist{i} = dealMultinomial([1-proba(1) proba(1)],getCellNb(model)) ;
        end
        conductivities{i,j} = distributeMicro(model,dist{i},microPatterns) ;
    end
    % QP approximation
    [KhomoQP,outQP{j}] = homogenize(pb,1,maxIter,stagnationTol,...
        conductivities(:,j)) ;
    % FE computations
    [KhomoFE,outFE{j}] = homogenize(pb,2,numel(outQP{j}.times),0,...
        outQP{j}.conductivities) ;
end
outQP = cat(1,outQP{:}) ;
outFE = cat(1,outFE{:}) ;

%% Post-processing

% matlab2tikz additional arguments
m2tikzArgin = {'showInfo',false,'showwarnings',false,'noSize',true,...
    'extraAxisOptions',{'ymajorgrids=true','grid style ={dotted}'}} ;

% Line specifications and legend
specP = {':','--','-'} ;
leg = cell(2,nParam) ;
spec = leg ;
for j = 1:nParam
    spec{1,j} = [specP{j},'xr'] ;
    spec{2,j} = [specP{j},'ob'] ;
    leg{1,j} = sprintf('FEM (contrast=%g)',contrastK(j)) ;
    leg{2,j} = sprintf('LR (contrast=%g)',contrastK(j)) ;
end

% Stagnation and time
figure
subplot(2,1,1)
hold on
for j = 1:nParam
    plot(2:numel(outFE(j).times),outFE(j).stagnations,spec{1,j})
    plot(2:numel(outQP(j).times),outQP(j).stagnations,spec{2,j})
end
hold off
ylabel('Stagnation')
legend(leg(:))
subplot(2,1,2)
hold on
for j = 1:nParam
    plot(1:numel(outFE(j).times),log10(outFE(j).times),spec{1,j})
    plot(1:numel(outQP(j).times),log10(outQP(j).times),spec{2,j})
end
hold off
ylabel('Time (log(s))')
legend(leg(:))
cleanfigure('pruneText',false); 
tikzName = sprintf('contrast%iD_stag-time.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})

% Average time with respect to contrast
figure
hold on
avgTimes = zeros(nParam,2) ;
for j = 1:nParam
    avgTimes(j,1) = mean(outFE(j).times) ;
    avgTimes(j,2) = mean(outQP(j).times) ;
end
plot(contrastK',log10(avgTimes(:,1)),'-xr')
plot(contrastK',log10(avgTimes(:,2)),'-ob')
hold off
ylabel('Average time (log(s))')
xlabel('Conductivity contrast')
legend({'FEM','LR'})
cleanfigure('pruneText',false); 
tikzName = sprintf('contrast%iD_avgTime.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})

% Homogenized coefficients convergence
figure
if caseChoice == 1
    subplot(2,1,1)
    hold on
    for j = 1:nParam
        feK11 = cellfun(@(c)c(1,1),outFE(j).effectiveConductivities) ;
        plot(1:numel(outFE(j).times),feK11/feK11(end),spec{1,j})
        qpK11 = cellfun(@(c)c(1,1),outQP(j).effectiveConductivities) ;
        plot(1:numel(outQP(j).times),qpK11/qpK11(end),spec{2,j})
    end
    hold off
    legend(leg(:))
    ylabel('\theta_{m}(K_{*11}^{400})/\theta_{20}(K_{*11}^{400})')
    subplot(2,1,2)
end
hold on
for j = 1:nParam
    feK22 = cellfun(@(c)c(2,2),outFE(j).effectiveConductivities) ;
    plot(1:numel(outFE(j).times),feK22/feK22(end),spec{1,j})
    qpK22 = cellfun(@(c)c(2,2),outQP(j).effectiveConductivities) ;
    plot(1:numel(outQP(j).times),qpK22/qpK22(end),spec{2,j})
end
hold off
legend(leg(:))
ylabel('\theta_{m}(K_{*22}^{400})/\theta_{20}(K_{*22}^{400})')
cleanfigure('pruneText',false); 
tikzName = sprintf('contrast%iD_K.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})