% Parameters

sampleSize = [20 20] ; % domain size for every sample
sampleNb = 30 ; % number of simulation per probability value
probaValues = [0.001 .01 0.05 .1:.1:.9 0.95 0.99 0.999] ;

%% Pre-processing

probaNb = numel(probaValues) ;
ranks = zeros(probaNb,sampleNb) ;
times = ranks ;

model = QPModel('order',2,...
    'cellNum',sampleSize,...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',false) ;

patterns = struct('name',{'uniform','disc','rectangle'},...
    'value',{1 99 99},...
    'size',{[] .25 sqrt([.5 .5])} ,...
    'center',{[] [.5 .5] [.5 .5]},...
    'offset',{[] [] []}) ;
patternsTable = [1 1 ; 0 0 ; 1 0] ;

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'probability',[0.9 0.1]) ;
KAss = assemble(KAss) ; % to get conductivity and its bounds

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'useCompression',true) ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

%% Processing

for n = 1:probaNb
    fprintf('%i_ Defect probability set to %.3f\n',n,probaValues(n));
    for m = 1:sampleNb
        fprintf(' Sample %i - ',m);
        dist = dealMultinomial([1-probaValues(n) probaValues(n)],...
            prod(sampleSize)) ;
       [pb,times(n,m)] = updateDistribution(pb,dist) ;
       [~,output] = solve(pb) ;
       times(n,m) = times(n,m) + output.time + output.initTime ;
       ranks(n,m) = output.iter ;
       fprintf('Rank %i - %s - error %.3g\n',ranks(n,m),...
           formatDuration(times(n,m)),output.error) ;
    end
end

%% Post-processing

avgRanks = mean(ranks,2) ;
avgTimes = mean(times,2) ;

% matlab2tikz additional arguments
m2tikzArgin = {'showInfo',false,'showwarnings',false,'noSize',true,...
    'extraAxisOptions',{'ymajorgrids=true','grid style ={dotted}'}} ;

figure
%subplot(2,2,1)
plot(probaValues(:)',avgRanks,'o-b')
xlabel('Defect probability')
ylabel('Average approximation rank')
title(sprintf('%i samples over %ix%i cells',sampleNb,sampleSize(1),...
    sampleSize(2)))

cleanfigure('pruneText',false); 
matlab2tikz('rank-wrt-proba_avg.tikz',m2tikzArgin{:})

%subplot(2,2,2)
% plot(probaValues(:)',avgTimes,'d-k')
% xlabel('Defect probability')
% ylabel('Average computation time')
% title(sprintf('%i samples over %ix%i cells',sampleNb,sampleSize(1),...
%     sampleSize(2)))
% 
figure% subplot(2,2,3)
plot(probaValues(:)',var(ranks,0,2),'o-b')
xlabel('Defect probability')
ylabel('Variance on average rank')
title(sprintf('%i samples over %ix%i cells',sampleNb,sampleSize(1),...
    sampleSize(2)))

cleanfigure('pruneText',false); 
matlab2tikz('rank-wrt-proba_var.tikz',m2tikzArgin{:})
% 
% subplot(2,2,4)
% plot(probaValues(:)',var(times,0,2),'d-k')
% xlabel('Defect probability')
% ylabel('Variance on average time')
% title(sprintf('%i samples over %ix%i cells',sampleNb,sampleSize(1),...
%     sampleSize(2)))