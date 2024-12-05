function [solution,output] = solveFE(pb,useMultiplier,meshMethod)
% [solution,output] = solveFE(pb,useMultiplier,meshMethod)

if nargin < 3
    meshMethod = 1 ;
    if nargin < 2
        useMultiplier = false ;
    end
end

feModel = [] ;
if isa(meshMethod,'MODEL')
    feModel = meshMethod ;
end

%% Parameters

qpModel = pb.model ;
cellNum = getCellNum(qpModel) ;
cellSize = getCellSize(qpModel) ;
cellNb = getCellNb(qpModel) ;
bc = pb.bc ;
cstNull = 'full' ;
output = QPProblem.feOutputStructure(1) ;
output.operators = struct('diffOp',{[]},'lhsBCOp',{[]},'srcOp', ...
    {[]},'rhsBCOp',{[]}, 'LagrangeOp',{[]}) ;

% Check multiplier used consistently, if used.
assert(~useMultiplier || (numel(bc)==1 && getType(bc(1))==1),...
    'Boundary conditions insconsistent with multiplier use.') ;

%% Model

% Mesh creation (choose method)
if isempty(feModel)
    [feModel,mesherTime] = buildFEModel(pb,meshMethod) ;
    feModel = concatgroupelem(feModel,[],'force'); % this should not be necessary (yet it is)
else
    mesherTime = 0 ;
end
periodicity = [isPeriodic(bc,1) isPeriodic(bc,2)] ;

if all(periodicity) && strcmp(cstNull,'1point')
    % Identify free corner node
    cornerCoord = [0 0 ; 1 0 ; 1 1 ; 0 1]*diag(cellSize.*cellNum) ;
    [~,cornerNodes] = ismember(cornerCoord,getcoord(getnode(feModel))) ;
    bcond = getbcond(feModel) ;
    free = ismember(cornerNodes,getddlfree(bcond)) ;
    % Set this node value to zero
    feModel = addcl(feModel,cornerNodes(free));
end


ifprint(pb.verboseGreedy,sprintf('%i periods meshed in %.2f s\n',cellNb,mesherTime))

%% Assembling

assemblerClock = tic ;

% Diffusion operator
coord = getcoord(getnode(feModel)) ;
relCoord = coord*diag(max(coord).^-1) ; % relative coordinates, in [0,1]^d
conductivityField = getCompleteK(pb,feModel) ;
a = BILINFORM(1,1,conductivityField,0) ;
a = setfree(a,0) ;
A = calc_matrix(a,feModel) ;
output.operators.diffOp = A ;

% Source operator
source = getCompleteSource(pb,feModel) ;
if ischar(source)
    if strcmp(source,'corrector1')
        b = -A*coord(:,1) ;
    elseif strcmp(source,'corrector2')
        b = -A*coord(:,2) ;
    elseif strcmp(source,'correctors')
        b = -A*coord ;
    else
        error('Unknown source')
    end
elseif isnumeric(source)
    b = calc_vector(setfree(LINFORM(0,source,0),0),feModel) ;
else
    error('Unknown source')
end
output.operators.srcOp = b ;

% Boundary conditions
if periodicity(1)
    domain = DOMAIN(2,[0 0],cellSize.*cellNum) ;
    feModel = addclperiodic(feModel,getedge(domain,1),getedge(domain,3),'u');
end
if periodicity(2)
    domain = DOMAIN(2,[0 0],cellSize.*cellNum) ;
    feModel = addclperiodic(feModel,getedge(domain,2),getedge(domain,4),'u');
end
output.operators.lhsBCOp = cell(numel(bc),1) ;
output.operators.rhsBCOp = cell(numel(bc),1) ;
for i = 1:numel(bc)
    if bc(i).type~=4 
        M = edgeOp(feModel) ;
        M = M{1} + M{2} + M{3} + M{4} ;
        bcLoc = bc(i).loc(relCoord) ;
        M(~bcLoc,:) = 0 ;
    end
    switch bc(i).type
        case 1
            if useMultiplier
                Ll = M ; % Lagrange LHS
                Lr = Ll*bc(i).eval(relCoord) ; % Lagrange RHS
                output.operators.lhsBCOp{i} = Ll ;
                output.operators.rhsBCOp{i} = Lr ;
            else
                N = edgeGradOp(feModel,conductivityField) ;
                N = N{1} + N{2} + N{3} + N{4} ;
                N(~bcLoc,:) = 0 ;
                penalty = bc(i).factor ;
                if isempty(penalty)
                    penalty = 1 ;
                end
                lhs = penalty*M - N - N' ;
                rhs = (penalty*M - N')*bc(i).eval(relCoord) ;
                A = A + lhs ;
                b = b + repmat(rhs,1,size(b,2)) ; % repmat if multiple src
                output.operators.lhsBCOp{i} = lhs ;
                output.operators.rhsBCOp{i} = rhs ;
            end
            
        case 2 % Neumann
            M = M*bc(i).eval(relCoord) ;
            output.operators.rhsBCOp{i} = M ;
            b = b + repmat(M,1,size(b,2))  ;
            
        case 3 % Robin
            coef = bc(i).factor ;
            if isnumeric(coef) && isscalar(coef)
                Mc = coef*M ;
            else
                Mc = edgeOp(feModel,coef) ;
                Mc = Mc{1} + Mc{2} + Mc{3} + Mc{4} ;
                Mc(~bcLoc,:) = 0 ;
            end
            M = M*bc(i).eval(relCoord) ;
            output.operators.lhsBCOp{i} = Mc ;
            output.operators.rhsBCOp{i} = M ;
            A = A + Mc ;
            b = b + repmat(M,size(b,2)) ;
            
        case 4 % Periodic
            % Do nothing: already strongly enforced.
    end
end
% Constant nullification
periodicOrNeumanBC = ismember([bc(:).type],[2 4]) ;
if all(periodicOrNeumanBC) % only Neumann and Periodic BC
    i = find(periodicOrNeumanBC,1) ; % first such BC
    output.operators.lhsBCOp{i} = calc_matrix(...
        setfree(BILINFORM(0,0,1),0),feModel) ;
    A = A + output.operators.lhsBCOp{i} ;
end
% Set empty cells to zero
emptyLHSBC = cellfun(@isempty,output.operators.lhsBCOp) ;
if any(emptyLHSBC)
    output.operators.lhsBCOp(emptyLHSBC) = ...
        {sparse(getnbddl(feModel),getnbddl(feModel))} ;
end
emptyRHSBC = cellfun(@isempty,output.operators.rhsBCOp) ;
if any(emptyRHSBC)
    output.operators.rhsBCOp(emptyRHSBC) = {sparse(getnbddl(feModel),1)} ;
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
ifprint(pb.verboseGreedy,sprintf('Finite element assembling in %.2f s\n',assemblerTime))

%% Resolution

solverClock = tic ;
solution = A\b ;
solverTime = toc(solverClock) ;
ifprint(pb.verboseGreedy,sprintf('Finite element computation in %.2f s\n',solverTime))

%% Post-Processing

if useMultiplier % Extract relevant operators and vectors
    range = 1:(size(A,1)-numel(Pj)) ;
    A = A(range,range) ;
    b = b(range,:) ;
    output.multiplier = sparse(size(b,1),size(b,2)) ;
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