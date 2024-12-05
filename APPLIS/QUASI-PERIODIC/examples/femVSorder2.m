cellNums = repmat([5 10 14 17 20 22 24 26 28 30 32]',1,2) ;
bcTypes = [5];% 1 2 3] ;
bcNames = {'P'};%,'D','N','R'} ;
bcSpecs = {'g','b','r','m'} ;
src = {'corrector1' ; 2.7 ; [] } ;
srcNames = {'Corrector','Uniform','Peak'} ;
srcSpecs = {'-d','--o',':^'} ;
N = size(cellNums,1) ;

model = QPModel('order',2,...
    'cellNum',cellNums(1,:),...
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

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[0.9 0.1]) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1') ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

M = numel(src)*numel(bcTypes) ; % shortcut
times = zeros(N,1+M) ;
ranks = zeros(N,M) ;
linespecs = cell(1,1+M) ;
linespecs{1} = '-dk' ;
legends = cell(1,1+M) ;
legends{1} = 'FEM-PBC-Corrector' ;
for icn = 1:N
    model = setCellNum(model,cellNums(icn,:)) ;
    dsz = getDomainSize(model) ;
    src{3} = @(x) exp(-10*sqrt((x(:,1)-dsz(1)/2).^2+(x(:,2)-dsz(2)/2).^2)) ;
    c = 2 ;
    for ibc = 1%:numel(bcTypes)
        for isrc = 1%:numel(src)
            bc = {QPBC('type',bcTypes(ibc),'model',model,'value',1.3,...
                'penalty',3.1)} ;
            KAss = updateModel(KAss,model) ;
            diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
                'BC',bc,...
                'source',src{isrc}) ;
            pb = QPDiffusionProblem('operatorAssembler',diffAss,...
                'tolerance',1e-3);
            [~,outQP,pb] = solve(pb,1:getOrder(pb)) ;
            if bcTypes(ibc)==5 && ischar(src{isrc})
                pbFE = pb ;
            end
            times(icn,c) = outQP.time + outQP.initTime;
            ranks(icn,c-1) = outQP.iter ;
            if icn == 1
                linespecs{c} = [srcSpecs{isrc},bcSpecs{ibc}] ;
                legends{c} = ['LR-',bcNames{ibc},'BC-',srcNames{isrc}] ;
            end
            c = c +1 ;
        end
    end
    [~,outFE] = solveFEM(pbFE,false,2) ;
    times(icn,1) = sum(outFE.time) ;
end

%%

% times(1:11,1) = [3.47; 30.12; 126.24; 263.27; 541.3276; 773.0571; ...
%     1068.10875; 1520.14; 2001.115; 2602.875; 3299.5] ;
cellNbs = prod(cellNums,2) ;

figure ; hold on
for i = 1:2
    plot(cellNbs,times(:,i),linespecs{i})
end
legend(legends(1:2))
xlabel('Number of cells')
ylabel('Time (s)')
hold off

figure ; hold on
for i = 2:1+M
    plot(cellNbs,times(:,i),linespecs{i})
end
legend(legends(2:end))
xlabel('Number of cells')
ylabel('Time (s)')
hold off

figure ; hold on
for i = 1:M
    plot(cellNbs,ranks(:,i),linespecs{i+1})
end
legend(legends(2:end))
xlabel('Number of cells')
ylabel('Rank')
hold off