function [Assembler,assemblerTime] = assembleMicroOperators(Assembler)
% [Assembler,assemblerTime] = assembleMicroOperators(Assembler)

%% Preparation

assemblerClock = tic ;

% For speed and clarity
cellModel = getCellModel(Assembler) ;
dim = getdim(cellModel) ; % must be 2
patternNb = getConductivityPatternNb(Assembler) ;
Ks = getConductivity(Assembler) ;
Ks = mat2cell(Ks.space.spaces{end},Ks.space.sz(end),...
    ones(1,Ks.space.dim(end)))' ; %TODO: this neglects the core

% Toward edges models and edge translation matrices
% Get edges
cellBoundary = create_boundary(cellModel,'withparent') ;
[cellEdges,normales] = getedges(getCellDomain(Assembler)) ;
%refNormales = [1 0 ; -1 0 ; 0 1 ; 0 -1] ; % for order reference
refNormales = [1 0 ; 0 1 ; -1 0 ; 0 -1] ; %TODO: test difference wrt orders

% Re-order Normales and Edges
[~,orderInd]=ismember(refNormales,cell2mat(normales)','rows') ;
normales = normales(orderInd) ;
cellEdges = cellEdges(orderInd) ;
edgeNb = numel(cellEdges) ; % must be 4 (=2*dim)

% Build edges models
cellEdgeModel = cell(1,edgeNb) ;
for e=1:edgeNb
    cellEdgeModel{e} = intersect(cellBoundary,cellEdges{e}) ;
end

% Edge to edge translation matrices
P{1} = calc_P_edges(cellModel,cellEdges{1},cellEdges{3}) ;
P{2} = calc_P_edges(cellModel,cellEdges{2},cellEdges{4}) ;

%% Diffusion

diffusionOperators = cell(patternNb,1) ;
for r = 1:patternNb
    diffusionForm = BILINFORM(1,1,Ks{r},0) ;
    diffusionForm = setfree(diffusionForm,0) ; % Unfree to apply no BC
    diffusionOperators{r} = calc_matrix(diffusionForm,cellModel) ;
end

%% Consistency

edgeGradientForm = BILINFORMBOUNDARY(0,1,1,0) ;
edgeGradientForm = setfree(edgeGradientForm,0) ; % Unfree to apply no BC
edgeGradientOp = cell(4,dim,patternNb) ;
for r = 1:patternNb % For each pattern
    for d = 1:dim % Per dim (in order) : N_0^ed N_0^(-ed) N_1^ed N_1^(-ed)
        NK = {normales{d}(1)*Ks{r} ; normales{d}(2)*Ks{r} } ;
        edgeGradientForm = setk(edgeGradientForm,NK) ;
        edgeGradientOp{1,d,r} = calc_matrix(edgeGradientForm,...
            cellEdgeModel{d},'parent',cellModel) ;
        edgeGradientOp{2,d,r} = -calc_matrix(edgeGradientForm,...
            cellEdgeModel{d+2},'parent',cellModel) ;
        edgeGradientOp{3,d,r} = P{d}*edgeGradientOp{1,d,r} ;
        edgeGradientOp{4,d,r} = P{d}'*edgeGradientOp{2,d,r} ;
    end
end

% Storage
% 4*dim*patternNb operators stored in patternNb columns.
% All operators consecutively, for every pattern K^Y_n. Operators in order:
% N_0^e1 N_0^-e1 N_1^e1 N_1^-e1 N_0^e2 N_0^-e2 N_1^e2 N_1^-e2
consistencyOperators = reshape(edgeGradientOp,[4*dim patternNb]) ;

%% Stabilisation

edgeForm = BILINFORMBOUNDARY(0,0,1) ;
edgeForm = setfree(edgeForm,0) ; % Unfree to apply no BC
edgeOp = cell(4,dim) ;
for d=1:dim
    edgeOp{1,d} = calc_matrix(edgeForm,cellEdgeModel{d},'parent',cellModel) ; % M_0^ed
    edgeOp{2,d} = calc_matrix(edgeForm,cellEdgeModel{d+2},'parent',cellModel) ; % M_0^{-ed}
    edgeOp{3,d} = P{d}*edgeOp{1,d} ; % M_1^{ed}
    edgeOp{4,d} = edgeOp{3,d}' ; % M_1^{-ed} (= {}^tM_1^{ed})
end

% Storage
% 4*dim operators stored in one column, in order:
% M_0^e1 M_0^{-e1} M_1^e1 M_1^{-e1} M_0^e2 M_0^-e2 M_1^e2 M_1^-e2
stabilisationOperators = edgeOp(:) ;

%% Constant nullification

switch getConstantNullification(Assembler) ;
    case '1point' % Lower Left node (supposedly)
        cstNullOperator = {sparse(1,1,1,getnbddl(cellModel),getnbddl(cellModel))} ;
    case 'full'
        cstNullOperator = {calc_matrix(setfree(BILINFORM(0,0,1),0),cellModel)} ;
    otherwise
        error('QPDiffusionAssembler: unknown constant nullification method')
end

%% RHS

source = getSource(Assembler) ;
calcRHSOperators = ~ischar(source) && isa(source.space,'TSpaceVectors') ;
if calcRHSOperators
    srcMicroVal = source.space.spaces{end} ;
    sourcePatternNb = size(srcMicroVal,2) ;
    rhsOperators = zeros(getnbddl(cellModel),sourcePatternNb) ;
    rhsLinForm = LINFORM(0,1,0) ;
    rhsLinForm = setfree(rhsLinForm,0) ; % Unfree to apply no BC
    for r=1:sourcePatternNb
        rhsLinForm = setk(rhsLinForm,srcMicroVal(:,r)) ;
        rhsOperators(:,r) = calc_vector(rhsLinForm,cellModel) ;
    end
else
    rhsOperators = [] ;
end

%% Operators storage

% Store as structure
microOperators = struct('diffusion',{diffusionOperators},'consistency',...
    {consistencyOperators},'stabilisation',{stabilisationOperators},...
    'cstNull',{cstNullOperator},'rhs',{rhsOperators}) ;

% Pass to attribute
Assembler = setMicroOperators(Assembler,microOperators) ;

assemblerTime = toc(assemblerClock) ;
end