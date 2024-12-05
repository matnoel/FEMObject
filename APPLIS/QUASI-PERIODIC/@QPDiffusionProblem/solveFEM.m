function [solution,output] = solveFEM(pb,useMultiplier,meshMethod)
% [solution,output] = solveFEM(pb,useMultiplier,meshMethod)

if nargin < 3
    meshMethod = 1 ;
    if nargin < 2
        useMultiplier = false ;
    end
end

%% Parameters

qpModel = getModel(pb) ;
cellNum = getCellNum(qpModel) ;
cellSize = getCellSize(qpModel) ;
cellNb = getCellNb(qpModel) ;
bc = getBC(pb) ;
source = getSource(pb) ;
cstNull = getConstantNullification(pb) ;
output = QPDiffusionProblem.feOutputStructure(1) ;
output.operators = struct('diffOp',{[]},'lhsBCOp',{[]},'srcOp', ...
    {[]},'rhsBCOp',{[]}, 'LagrangeOp',{[]}) ;

% Check multiplier used consistently, if used.
assert(~useMultiplier || (numel(bc)==1 && getType(bc{1})==1),...
    'Boundary conditions insconsistent with multiplier use.') ;

%% Model

% Mesh creation (choose method)
[feModel,mesherTime] = buildFEModel(qpModel,meshMethod) ;

if isPeriodic(pb) && strcmp(cstNull,'1point')
    % Identify free corner node
    cornerCoord = [0 0 ; 1 0 ; 1 1 ; 0 1]*diag(cellSize.*cellNum) ;
    [~,cornerNodes] = ismember(cornerCoord,getcoord(getnode(feModel))) ;
    bc = getbcond(feModel) ;
    free = ismember(cornerNodes,getddlfree(bc)) ;
    % Set this node value to zero
    feModel = addcl(feModel,cornerNodes(free));
end


ifprint(pb,sprintf('%i periods meshed in %.2f s\n',cellNb,mesherTime))

%% Assembling

assemblerClock = tic ;

% Diffusion operator
if isempty(getPatches(pb))
    conductivityField = untensorize(qpModel,getConductivity(pb)) ;
else % get real conductivity
    conductivityField = untensorize(qpModel,getConductivity(getPatches(pb))) ;
end
coord = getcoord(getnode(feModel)) ;
qpCoord = smooth(qpModel,getDomainCoord(qpModel)) ;
conductivityField = relativeSort(conductivityField,qpCoord,coord) ;
a = BILINFORM(1,1,conductivityField,0) ;
a = setfree(a,0) ;
A = calc_matrix(a,feModel) ;
output.operators.diffOp = A ;

% Source operator
switch class(source)
    case 'char'
        if strcmp(source,'corrector1')
            b = A*coord(:,1) ;
        elseif strcmp(source,'corrector2')
            b = A*coord(:,2) ;
        elseif strcmp(source,'correctors')
            b = A*coord ;
        else
            error('QPDiffusionPoblem.solveFEM: unknown source')
        end
    case 'double'
        b = relativeSort(source,qpCoord,coord) ;
        b = calc_vector(setfree(LINFORM(0,b,0),0),feModel) ;
    case 'TuckerLikeTensor'
        b = untensorize(qpModel,source) ;
        b = relativeSort(b,qpCoord,coord) ;
        b = calc_vector(setfree(LINFORM(0,b,0),0),feModel) ;
    otherwise
        error('QPDiffusionPoblem.solveFEM: unknown source')
end
output.operators.srcOp = b ;

% Boundary conditions
if isPeriodic(pb,1)
    domain = DOMAIN(2,[0 0],cellSize.*cellNum) ;
    feModel = addclperiodic(feModel,getedge(domain,1),getedge(domain,3),'u');
end
if isPeriodic(pb,2)
    domain = DOMAIN(2,[0 0],cellSize.*cellNum) ;
    feModel = addclperiodic(feModel,getedge(domain,2),getedge(domain,4),'u');
end
bc = getBC(pb) ;
output.operators.lhsBCOp = cell(numel(bc),1) ;
output.operators.rhsBCOp = cell(numel(bc),1) ;
for i = 1:numel(bc)
    val = getValue(bc{i}) ;
    if iscell(val) % Cauchy
        val = val{2} ;
    end
    val = untensorize(qpModel,val) ;
    val = relativeSort(val,qpCoord,coord) ;
    switch getType(bc{i})
        case 1
            if useMultiplier
                Ll = assembleEdgeOperator(bc{i},feModel) ; % Lagrange LHS
                Lr = Ll*val ; % Lagrange RHS
                output.operators.lhsBCOp{i} = Ll ;
                output.operators.rhsBCOp{i} = Lr ;
            else
                [lhs,rhs] = assembleNitsche(bc{i},feModel,conductivityField) ;
                A = A + lhs ;
                b = b + repmat(rhs,1,size(b,2)) ; % repmat if multiple src
                output.operators.lhsBCOp{i} = lhs ;
                output.operators.rhsBCOp{i} = rhs ;
            end
            
        case 2 % Neumann
            output.operators.rhsBCOp{i} = assembleEdgeOperator(bc{i},...
                feModel)*val ;
            b = b + repmat(output.operators.rhsBCOp{i},1,size(b,2))  ;
            
        case 3 % Robin
            coef = untensorize(qpModel,getPenalty(bc{i})) ;
            coef = relativeSort(coef,qpCoord,coord) ;
            output.operators.lhsBCOp{i} = assembleEdgeOperator(bc{i},...
                feModel,full(coef)) ;
            output.operators.rhsBCOp{i} = assembleEdgeOperator(bc{i},...
                feModel)*val ;
            A = A + output.operators.lhsBCOp{i} ;
            b = b + repmat(output.operators.rhsBCOp{i},size(b,2)) ;
            
        case 4 % Cauchy
            % Dirichlet
            [lhs,rhs] = assembleNitsche(bc{i},feModel,conductivityField) ;
            % Neumann
            rhs = rhs + assembleEdgeOperator(bc{i},feModel)*val ;
            % Add and store
            A = A + lhs ;
            b = b + repmat(rhs,1,size(b,2)) ;
            output.operators.lhsBCOp{i} = lhs ;
            output.operators.rhsBCOp{i} = rhs ;
            
        case 5 % Periodic
            % Do nothing: already strongly enforced.
    end
end
% Constant nullification
NeumannOrPeriodicBC = ismember(cellfun(@getType,bc),[2 5]) ;
if all(NeumannOrPeriodicBC) % only Neumann and Periodic BC
    i = find(NeumannOrPeriodicBC,1) ; % first such BC
    if strcmp(cstNull,'full')
        output.operators.lhsBCOp{i} = calc_matrix(...
            setfree(BILINFORM(0,0,1),0),feModel) ;
        A = A + output.operators.lhsBCOp{i} ;
    elseif strcmp(cstNull,'1point')
        A(1,1) = A(1,1) + 1 ;
        output.operators.lhsBCOp{i} = sparse(1,1,1,size(A,1),size(A,2)) ;
    end
end
% Set empty cells to zero
emptyLHSBC = cellfun(@isempty,output.operators.lhsBCOp) ;
if any(emptyLHSBC)
    output.operators.lhsBCOp(emptyLHSBC) = {zeros(getnbddl(feModel))} ;
end
emptyRHSBC = cellfun(@isempty,output.operators.rhsBCOp) ;
if any(emptyRHSBC)
    output.operators.rhsBCOp(emptyRHSBC) = {zeros(getnbddl(feModel),1)} ;
end

A = freematrix(feModel,A) ;
b = freevector(feModel,b) ;
if useMultiplier
    [Pi,Pj] = find(Ll) ;
    Pj = unique(Pj) ;
    A = [A -Ll(:,Pj) ; Ll(unique(Pi),:) zeros(numel(Pj))] ;
    b = [b ; repmat(Lr(Pj),1,size(b,2))] ;
end

assemblerTime = toc(assemblerClock) ;
ifprint(pb,sprintf('Finite element assembling in %.2f s\n',assemblerTime))

%% Resolution

solverClock = tic ;
solution = A\b ;
solverTime = toc(solverClock) ;
ifprint(pb,sprintf('Finite element computation in %.2f s\n',solverTime))

%% Post-Processing

if useMultiplier % Extract relevant operators and vectors
    range = 1:(size(A,1)-numel(Pj)) ;
    A = A(range,range) ;
    b = b(range,:) ;
    output.multiplier = zeros(size(b)) ;
    output.multiplier(Pj,:) = solution(setdiff(1:numel(solution),range),:) ;
    solution = solution(range,:) ;
end

% Residual error
res = A*solution-b ;
output.error = diag(res'*res)./diag(b'*b) ; % no 'norm' to handle
% multiple src case
output.residual = unfreevector(feModel,res) ;

% Storage
output.time = [mesherTime assemblerTime solverTime] ;
output.model = feModel ;
output.lhs = unfreematrix(feModel,A) ;
output.rhs = unfreevector(feModel,b) ;
solution = unfreevector(feModel,solution) ;

end

function [lhs,rhs] = assembleNitsche(bc,feModel,K)

qpModel = getModel(bc) ;
cellNum = getCellNum(qpModel) ;
coord = getcoord(getnode(feModel)) ;
dSz = max(coord) ;
domain = DOMAIN(2,[0 0],dSz) ;
mesoCoord = formatIndex(2,cellNum,mesoMicroCoord(qpModel,coord)) ;

% Boundary
boundary = create_boundary(feModel,'withparent') ;
[edges,normales] = getedges(domain) ;
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
normales = normales(orderInd) ;
edges = edges(orderInd) ;
bCells = getCells(bc) ;

% Edge operators
boundaryForm = setfree(BILINFORMBOUNDARY(0,0,1,0),0) ;
edgeOp = zeros(size(coord,1)) ;
boundaryGradForm = setfree(BILINFORMBOUNDARY(0,1,1,0),0) ;
edgeGradOp = zeros(size(coord,1)) ;
for e = 1:4
    % Select relevant DoFs
    adir = 1+mod(e-1,2) ; % absolute direction (1 or 2)
    sdir = sign(2.5-e) ; % direction sign (-1 or 1)
    edgeCoord = max(0,sdir)*dSz(adir) ;
    loc = coord(:,adir) == edgeCoord ; % DoFs on relevant edge
    loc = loc & ismember(mesoCoord,formatIndex(2,cellNum,bCells{e})) ; % restrict to bc cells
    % Localize then assemble edge operators
    edgeModel = intersect(boundary,edges{e}) ;
    edgeForm = setk(boundaryForm,double(loc)) ;
    edgeOp = edgeOp + calc_matrix(edgeForm,edgeModel,'parent',feModel) ;
    Kloc = K.*loc ;
    NK = {normales{e}(1)*Kloc ; normales{e}(2)*Kloc} ;
    edgeGradForm = setk(boundaryGradForm,NK) ;
    edgeGradOp = edgeGradOp + calc_matrix(edgeGradForm,edgeModel,'parent',feModel) ;
end

% Assemble
penalty = getPenalty(bc) ;
val = getValue(bc) ;
if iscell(val) % To process Cauchy QPBC
    val = val{1} ;
end
val = untensorize(qpModel,val) ;
qpCoord = smooth(qpModel,getDomainCoord(qpModel)) ;
val = relativeSort(val,qpCoord,coord) ;
lhs = penalty*edgeOp - edgeGradOp - edgeGradOp' ;
rhs = (penalty*edgeOp - edgeGradOp')*val ;

end

function lhs = assembleEdgeOperator(bc,feModel,coef)

% This is an excerpt from assembleNitsche, to assemble only "edgeOp".
% We leave assembleNitsche as a whole for performance; breaking it down
% into smaller functions would require to repeat some operations.

if nargin < 3
    coef = ones(getnbddl(feModel),1) ;
end

qpModel = getModel(bc) ;
cellNum = getCellNum(qpModel);
coord = getcoord(getnode(feModel)) ;
dSz = max(coord) ;
domain = DOMAIN(2,[0 0],dSz) ;
mesoCoord = formatIndex(2,cellNum,mesoMicroCoord(qpModel,coord)) ;

% Boundary model, sort edges
boundary = create_boundary(feModel,'withparent') ;
[edges,normales] = getedges(domain) ;
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ;
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
edges = edges(orderInd) ;
bCells = getCells(bc) ;

% Edge operator
boundaryForm = setfree(BILINFORMBOUNDARY(0,0,1,0),0) ;
edgeOp = zeros(size(coord,1)) ;
for e = 1:4
    % Select relevant DoFs
    adir = 1+mod(e-1,2) ; % absolute direction (1 or 2)
    sdir = sign(2.5-e) ; % direction sign (-1 or 1)
    edgeCoord = max(0,sdir)*dSz(adir) ;
    loc = coord(:,adir) == edgeCoord ; % DoFs on relevant edge
    loc = loc & ismember(mesoCoord,formatIndex(2,cellNum,bCells{e})) ; % restrict to bc cells
    
    edgeForm = setk(boundaryForm,coef.*loc) ;
    edgeModel = intersect(boundary,edges{e}) ;
    edgeOp = edgeOp + calc_matrix(edgeForm,edgeModel,'parent',feModel) ;
end

lhs = edgeOp ;

end