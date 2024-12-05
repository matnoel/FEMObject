function [sol,out,pb] = discBenchmark(order)

if nargin < 1
    order = 2 ;
end

cellNum = 10*[1 1] ;
model = QPModel('order',order,...
    'cellNum',cellNum,...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','disc'},...
    'value',{100 1},...
    'size',{[] .25} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
discCell = sub2ind(cellNum,ceil(cellNum(2)/2),ceil(cellNum(1))/2) ;
distribution = {(1:prod(cellNum))' discCell} ;
patternsTable = [1 1 ; 0 1] ;

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'distribution',distribution) ;

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'constantNullification','full') ;
diffAss = setPenalty(diffAss,50) ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

[sol,out,pb] = solve(pb,1:getOrder(pb)) ;

end