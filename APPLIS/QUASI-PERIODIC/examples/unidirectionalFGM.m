penalty = 1e3 ;
isPeriodic = true ;

model = QPModel('order',2,...
    'cellNum',[10 1],...
    'cellSize',[1 5],...
    'elementSize',[.05 .25],...
    'tolSVD',1e-6,...
    'verbose',true) ;

patterns = struct('name',{'uniform','bar'},...
    'value',{1 99},...
    'size',{[] [.2 1]} ,...
    'center',{[] [.5 .5]},...
    'offset',{[] []}) ;
patternsTable = [1 ; 1] ;

KAss = QPConductivityAssembler('model',model,...
    'patterns',patterns,...
    'patternsTable',patternsTable,...
    'distribution',{(1:getCellNb(model))'}) ;
KAss = assemble(KAss) ; % to get conductivity and its bounds

fgmPatterns = struct('name',{'bar','bar'},...
    'value',{1 1},...
    'size',{[.4 1] [.4 1]},...
    'center',{[.2 1] [.6 1]},...
    'offset',{[] []}) ;
fgmFields = drawCellPattern(getCellCoord(model),fgmPatterns) ;
J = TuckerLikeTensor.ones(tensorSize(model)) ;
J.space.spaces{end} = [J.space.spaces{end} fgmFields] ;
% J.core = FullTensor(ones(J.space.dim),getOrder(model),J.space.dim) ;
dilatation = [1 ; (1:getCellNb(model))'] ;
J.space.spaces{1}(:,2) = dilatation(1:end-1).*J.space.spaces{1}(:,1) ;
J.space.spaces{1}(:,3) = dilatation(2:end).*J.space.spaces{1}(:,1) ;
J.core = DiagonalTensor(ones(3,1),getOrder(model)) ;
J = updateAllProperties(J) ;
% KAss = setConductivity(KAss,getConductivity(KAss).*J) ;
K = getConductivity(KAss) ;
KJ = K.*J ;

% SIP Operator
cells = repmat(formatIndex(getOrder(model),getCellNum(model),...
    (1:getCellNb(model))'),1,2) ;
diffOp = bilinFormOperator(model,1,cells,KJ) ;
consOp = consistencyOperator(model,KJ,isPeriodic) ;
stabOp = penalty*stabilisationOperator(model,J,isPeriodic) ; % mesh measure missing
cstNull = bilinFormOperator(model,0) ;
sipOp = diffOp + stabOp - consOp - consOp' + cstNull ;

% Source Operator
src = getCoord(model) ;
src = (diffOp-consOp)*src{1} ; % corrector 1
srcOp = bilinFormOperator(model,0,cells,1)*src ; 

% Solver
greedy = GreedyLinearSolver(sipOp,srcOp,'tolerance',1e-2) ;
greedy.solutionUpdater = TuckerLikeSolutionUpdater(sipOp,srcOp,...
                'tolerance',1e-2,...
                'maxIterations',1,...
                'updateCore',false,...
                'updateTSpace',true,...
                'updateTSpaceDimensions',1:getOrder(model)) ;
            
% Resolution
[sol,out] = solve(greedy) ;

%% Post-processing

coord = getCoord(model) ;
x = untensorize(model,J.*coord{1},true) ;
y = untensorize(model,coord{2},true) ;
z = untensorize(model,sol,true) ;

figure ; plot(model,z,'coord',[x y])

