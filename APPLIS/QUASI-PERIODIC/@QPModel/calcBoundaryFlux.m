function flux = calcBoundaryFlux(model,v,microOp,bCells,K)
% flux = calcBoundaryFlux(model,v,microOp,bCells,K)
% Computes flux of v outward cells bCells, with respect to conductivity K.

if nargin < 5
    tsz = tensorSize(model) ;
    K = TuckerLikeTensor.ones(tsz,ones(1,numel(tsz))) ;
    if nargin < 4
        bCells = boundaryCells(model) ;
    end
end

% Preparations
if isnumeric(bCells)
        bCells = boundaryCells(model,bCells) ;
end
order = getOrder(model) ;
mesoOp = mesoIndicator(model,bCells,0) ;
microEdgeOp = microOp.stabilisation([1 5 2 6]') ;
microEdgeGradOp = assembleMicroGradOperator(model,K) ;
Ksp = K.space ; % for speed and brevity

%% Build operators
for d = 1:numel(bCells)
    % edgeOp along direction d
    % edgeOp space
    MSpace = cell(order,1) ;
    MSpace{end} = microEdgeOp(d) ;
    for o = 1:order-1
        MSpace{o} = mat2cell(mesoOp{o,d},size(mesoOp{o,d},1),...
            ones(1,size(mesoOp{o,d},2))) ;
        MSpace{o} = cellfun(@diag,MSpace{o}(:),'UniformOutput',0) ;
    end
    MSpace = TSpaceOperators(MSpace) ;
    
    % edgeOp core
    if order == 2 % each indicator is elementary
        MCore = DiagonalTensor(ones(MSpace.dim(1),1),v.order) ;
    else % order 3
        coreSz = MSpace.dim(1,:) ;
        indicatorRank = cellfun(@(x) size(x,2),mesoOp(1,:)) ;
        MCore = zeros(coreSz) ;
        for i = 1:coreSz(end)
            currentRank = 1+sum(indicatorRank(1:i-1)) ;
            range = currentRank:(currentRank+indicatorRank(i)-1) ;
            MCore(range,range,i) = eye(indicatorRank(i)) ;
        end
        MCore = FullTensor(MCore,order,MSpace.dim) ;
    end
    
    % edgeOp tensor
    M = TuckerLikeTensor(MCore,MSpace) ;
    
    % edgeGradOp along direction d
    NSpace = cell(order,1) ;
    NSpace{end} = microEdgeGradOp(d,:)' ;
    for o = 1:order-1
        NSpace{o} = cell(MSpace.dim(o),Ksp.dim(o)) ;
        for r = 1:Ksp.dim(o)
            for i = 1:MSpace.dim(o)
                NSpace{o}{i,r} = MSpace.spaces{o}{i}*diag(Ksp.spaces{o}(:,r)) ;
            end
        end
        NSpace{o} = NSpace{o}(:) ; % Sorted by K pattern last
    end
    NSpace = TSpaceOperators(NSpace) ;
    NCore = superkron(double(K.core),double(MCore)) ;
    NCore = FullTensor(NCore,order,MCore.sz.*K.core.sz) ;
    N = TuckerLikeTensor(NCore,NSpace) ;
    
    % Add
    if d == 1
        edgeOp = M ;
        edgeGradOp = N ;
    else
        edgeOp = edgeOp + M ;
        edgeGradOp = edgeGradOp + N ;
    end
end

%% Resolution

rhs = edgeGradOp*v ;
A = doubleQP(edgeOp) ;
b = doubleQP(rhs) ;
flux = zeros(size(b)) ;
[Pi,Pj] = find(A) ;
P = unique(Pj) ;
A = A(:,P) ;
A = A(unique(Pi),:) ;
b = b(P) ;
flux(P) = A\b ;
flux = tensorize(model,flux) ;

% [edgeOp,rhs] = convertTensors(edgeOp,rhs) ;
% tol = getTolSVD(model) ;
% bCellNb = sum(cellfun(@(x) size(x,1),bCells)) ;
% initialPoint = TuckerLikeTensor.ones(tensorSize(model)) ;
% %TODO: is it relevant to use tolSVD ? At least it is low.
% 
% localSolver = RankOneALSLinearSolver(edgeOp,rhs,'maxIterations',10,...
%     'stagnation',tol,'display',false,'x0',initialPoint) ;
% updater = TuckerLikeSolutionUpdater(edgeOp,rhs,'tolerance',tol, ...
%     'maxIterations',1,'stagnation',tol/10,'display',false,...
%     'updateCore',true,'updateTSpace',true,'updateDimTSpace',1:order-1) ;
% greedy = GreedyLinearSolver(edgeOp,rhs,'tolerance',tol,'stagnation',tol/10,...
%     'localSolver',localSolver,'solutionUpdater',updater,...
%     'maxIterations',bCellNb,'minimizeResidual',false,... % ?
%     'checkResidual',1,'display',false,'x0',initialPoint);
% 
% flux = solve(greedy) ;
end

function op = assembleMicroGradOperator(model,K)

microModel = getCellModel(model) ;
boundary = create_boundary(microModel,'withparent') ;
[edges,normales] = getedges(getCellDomain(model)) ;
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;

% Re-order Normales and Edges
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
normales = normales(orderInd) ;
edges = edges(orderInd) ;
edgeNb = numel(edges) ; % must be 4 (=2*dim)

% Build edges models
cellEdgeModel = cell(1,edgeNb) ;
for e=1:edgeNb
    cellEdgeModel{e} = intersect(boundary,edges{e}) ;
end

edgeGradientForm = BILINFORMBOUNDARY(0,1,1,0) ;
edgeGradientForm = setfree(edgeGradientForm,0) ; % Unfree to apply no BC
op = cell(edgeNb,K.space.dim(end)) ;
Ks = K.space.spaces{end} ;
dup = [] ;
for r = 1:K.space.dim(end) % For each pattern
    for e = 1:edgeNb % Per dim (in order) : N_0^ed N_0^(-ed)
        NK = {normales{e}(1)*Ks(:,r) ; normales{e}(2)*Ks(:,r) } ;
        edgeGradientForm = setk(edgeGradientForm,NK) ;
        op{e,r} = calc_matrix(edgeGradientForm,...
            cellEdgeModel{e},'parent',microModel) ;
        if e > 1
           dup = union(dup, intersect(find(op{e-1,r}),find(op{e,r})) ) ; 
        end
    end
    dup = union(dup, intersect(find(op{4,r}),find(op{1,r})) ) ; 
end

% % DEBUG
% Apply .5 factor to corner DoF counted twice.
halfactor = zeros(size(op{1})) ;
halfactor(dup) = .5 ;
% op = cellfun(@(x)x.*halfactor,op,'UniformOutput',false) ;
end