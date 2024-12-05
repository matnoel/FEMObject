sol = cell(6,1) ;
soltr = sol ;
dists = sol ;

model = QPModel('order',3,...
    'cellNum',[6 6],...
    'cellSize',[1 1],...
    'elementSize',[.05 .05],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','rectangle'},...
    'value',{1 99},...
    'size',{[] [.25 .25]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 1 ; 0 1] ;

cellList = (1:getCellNb(model))';
faulty = sub2ind(getCellNum(model),2,2) ;
dists{1} = {setdiff(cellList,faulty) faulty} ;
KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'distribution',dists{1}) ;
KAss = assemble(KAss) ; % to get conductivity and its bounds

diffAss = QPDiffusionAssembler('conductivityAssembler',KAss,...
    'BC','PBC',...
    'source','corrector1',...
    'useCompression',true) ;

pb = QPDiffusionProblem('operatorAssembler',diffAss,...
    'tolerance',1e-3);

dists{2}{2} = [faulty ; sub2ind(getCellNum(model),5,5)] ;
dists{2}{1} = setdiff(cellList,dists{2}{2}) ;
dists{3}{2} = [faulty ; faulty+6] ;
dists{3}{1} = setdiff(cellList,dists{3}{2}) ;
dists{4}{2} = [faulty ; faulty+6 ; sub2ind(getCellNum(model),5,5)] ;
dists{4}{1} = setdiff(cellList,dists{4}{2}) ;
dists{5}{2} = [faulty ; faulty+6 ; faulty+1] ;
dists{5}{1} = setdiff(cellList,dists{5}{2}) ;
dists{6}{2} = [faulty ; faulty+6 ; faulty+12] ;
dists{6}{1} = setdiff(cellList,dists{6}{2}) ;
[sol,out,pb] = multiSolve(pb,dists) ;

%%
tr=Truncator('tolerance',getTolSVD(model));

for i=1:numel(sol)
    soltr{i} = tr.truncate(sol{i});
    figure(i)
    plot(model,sol{i})
    title(sprintf('rank %i (%i)',sum([out(1:i).iter]),soltr{i}.space.dim(1)))
end