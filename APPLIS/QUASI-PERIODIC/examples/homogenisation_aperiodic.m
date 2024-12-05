% Parameters

caseChoice = 1 ; % 1 unidirectional | 2 bidirectional
maxIter = 20 ;
stagnationTol = 0 ;
approxTol = 1e-2 ;
mesoSize = 49 ;
contrastK = 100 ;
proba = .1 ;
valFun = @(n)rand(n) ; % @(n)1+.1*randn(n,1)
ptrnSzVal = [.2 .4 .6 .8] ;
ptrnFun = @(n) ptrnSzVal(ceil(rand(n,1)*numel(ptrnSzVal)))' ; % ensure column
% Change legend at the end if those last two are changed.

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
    szFun = @(sz) [sz 1] ;
    
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
    szFun = @(sz) [sz sz] ;
end

KAss = QPConductivityAssembler('model',model,...
                               'patterns',patterns,...
                               'patternsTable',patternsTable,...
                               'probability',[1-proba proba]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
                               'BC','PBC',...
                               'source','corrector1') ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,'tolerance',approxTol(1)) ;

randPtrn = [false true true] ;
randVal = [true false true] ;
nParam = numel(randPtrn) ;
cellNb = getCellNb(model) ;
cellCoord = getCellCoord(model) ;
conductivities = cell(maxIter,nParam) ;
tr = Truncator('tolerance',getTolSVD(model)) ;
patternsSizes = cellfun(@(x)ptrnFun(cellNb),cell(maxIter,1),...
    'UniformOutput',false) ;
randValues = cellfun(@(x)valFun(cellNb*[1 1]),cell(maxIter,1),...
    'UniformOutput',false) ;
outQP = cell(1,nParam) ;
outFE = cell(1,nParam) ;
for j = 1:nParam
    % Generate conductivities
    for i = 1:maxIter
        if randPtrn(j)
            ptrnSz = num2cell(patternsSizes{i}) ;
            ptrnSz = cellfun(szFun,ptrnSz,'UniformOutput',false) ;
            ptrns = patterns([1 ; 2*ones(cellNb,1)]) ;
            [ptrns(2:end).size] = deal(ptrnSz{:}) ;
            dist = [{(1:getCellNb(model))'} ; num2cell(1:cellNb)'] ;
        else
            ptrns = patterns ;
%             dist = dealMultinomial([proba 1-proba],cellNb) ;
            dist = repmat({(1:getCellNb(model))'},2,1) ;% no patternTable here
        end
        
        microPatterns = drawCellPattern(cellCoord,ptrns) ;
        K = distributeMicro(model,dist,microPatterns) ;
        
        if randVal(j)
            val = randValues{i}(:,1:K.space.dim(1)-1) ;
            K.space.spaces{1}(:,2:end) = K.space.spaces{1}(:,2:end).*val ;
        end
        
        K = tr.truncate(K) ;
        conductivities{i,j} = K ;
    end
    
    % QP approximation
    [KhomoQP,outQP{j}] = homogenize(pb,1,maxIter,stagnationTol,...
        conductivities(:,j)) ;
    % FE computations
%     [KhomoFE,outFE{j}] = homogenize(pb,2,numel(outQP{j}.times),0,...
%         outQP{j}.conductivities) ;
end
outQP = cat(1,outQP{:}) ;
outFE = cat(1,outFE{:}) ;

%% Post-processing

% matlab2tikz additional arguments
m2tikzArgin = {'showInfo',false,'showwarnings',false,'noSize',true,...
    'extraAxisOptions',{'ymajorgrids=true','grid style ={dotted}'}} ;

% Line specifications and legend
specVal = 'o' ; % 'o' if randVal 
specPtrn = '-' ; % '--' if randPtrn 
nameVal = {'values'} ;
namePtrn = {'patterns'} ;
leg = cell(2,nParam) ;
spec = leg ;
for j = 1:nParam
    spec{1,j} = ['-',specPtrn(randPtrn(j)),specVal(randVal(j)),'r'] ;
    spec{2,j} = ['-',specPtrn(randPtrn(j)),specVal(randVal(j)),'b'] ;
    leg{1,j} = sprintf('FEM (%s + %s)',nameVal{randVal(j)},...
        namePtrn{randPtrn(j)}) ;
    leg{2,j} = sprintf('LR (%s + %s)',nameVal{randVal(j)},...
        namePtrn{randPtrn(j)}) ;
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
tikzName = sprintf('aperiodic%iD_stag-time.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})

% Average times
% avgTimesAperiodic = zeros(nParam,2) ;
% for j = 1:nParam
%     avgTimesAperiodic(j,1) = mean(outQP(j).times) ;
%     avgTimesAperiodic(j,2) = mean(outFE(j).times) ;
% end
% [uPtrnNb,iu] = unique(ptrnNb) ;
% iu = ismember((1:numel(ptrnNb))',iu) ;
% leg2 = cell(4,1) ;
% figure
% hold on
% plot(uPtrnNb,avgTimes(iu,1),'--b')
% leg2{1} = sprintf('LR (%s)',randName{1}) ;
% plot(uPtrnNb,avgTimes(iu,2),'--r')
% leg2{2} = sprintf('FEM (%s)',randName{1}) ;
% plot(uPtrnNb,avgTimes(~iu,1),'-b')
% leg2{3} = sprintf('LR (%s)',randName{2}) ;
% plot(uPtrnNb,avgTimes(~iu,2),'-r')
% leg2{4} = sprintf('FEM (%s)',randName{2}) ;
% hold off
% ylabel('Average time (s)')
% xlabel('Pattern number')
% legend(leg2)

% Homogenized coefficients convergence (V1)
figure
if caseChoice == 1
    subplot(2,1,1)
    hold on
    for j = 1:nParam
        plot(1:numel(outFE(j).times),cellfun(@(c)c(1,1),...
            outFE(j).effectiveConductivities),spec{1,j})
        plot(1:numel(outQP(j).times),cellfun(@(c)c(1,1),...
            outQP(j).effectiveConductivities),spec{2,j})
    end
    hold off
    ylabel('K_{11}')
    legend(leg(:))
    subplot(2,1,2)
end
hold on
for j = 1:nParam
    plot(1:numel(outFE(j).times),cellfun(@(c)c(2,2),...
        outFE(j).effectiveConductivities),spec{1,j})
    plot(1:numel(outQP(j).times),cellfun(@(c)c(2,2),...
        outQP(j).effectiveConductivities),spec{2,j})
end
hold off
ylabel('K_{22}')
legend(leg(:))
cleanfigure('pruneText',false); 
tikzName = sprintf('aperiodic%iD_K.tikz',caseChoice) ;
matlab2tikz(tikzName,m2tikzArgin{:})

% % Homogenized coefficients convergence (V2)
% figure
% sf = subplot_format(2*nParam) ;
% for j = 1:nParam
%     subplot(sf(1),sf(2),j)
%     hold on
%     plot(1:numel(outFE(j).times),cellfun(@(c)c(1,1),...
%         outFE(j).effectiveConductivities),spec{1,j})
%     plot(1:numel(outQP(j).times),cellfun(@(c)c(1,1),...
%         outQP(j).effectiveConductivities),spec{2,j})
%     legend(leg(:,j))
%     ylabel('K_{11}')
%     hold off
%     
%     subplot(sf(1),sf(2),j+nParam)
%     hold on
%     plot(1:numel(outFE(j).times),cellfun(@(c)c(2,2),...
%         outFE(j).effectiveConductivities),spec{1,j})
%     plot(1:numel(outQP(j).times),cellfun(@(c)c(2,2),...
%         outQP(j).effectiveConductivities),spec{2,j})
%     legend(leg(:,j))
%     ylabel('K_{22}')
%     hold off
% end