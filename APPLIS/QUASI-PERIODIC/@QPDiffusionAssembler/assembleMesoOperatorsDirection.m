function mesoOperators = assembleMesoOperatorsDirection(assembler,dir)
% mesoOperators = assembleMesoOperatorsOrder(Assembler,dir)

order = getOrder(assembler) ;

% Cell connectivity matrix : cell i connected to cell i+1 (one-way only)
cellConn = calcMesoConnectivity(assembler,dir) ;
cellNb = size(cellConn,1) ; % total cell number for order 2, only along dir for order 3
% Diagonal matrices of cell connectivity rows and columns sum

%% Stabilisation

cellSize = getCellSize(assembler) ;
orthDir = 1+mod(dir+2,2) ; % 2 if direction=1, 1 if direction=2

if getUseStabilisationWeights(assembler)
    stabWeights = getStabilisationWeights(assembler,min(dir,order-1)) ;
    if ~iscell(stabWeights)
        stabWeights = {stabWeights} ;
    elseif size(stabWeights,1)~=1 % ensure single row format
        stabWeights = stabWeights(:)' ;
    end
else
    stabWeights = {ones(size(cellConn))} ;
end

% For stabilisation operator, apply weight scaling to connectivity matrix
IdW = cellfun(@(x) eye(size(x)).*x,stabWeights,...
    'UniformOutput',false) ;
cellConnW = cellfun(@(x) cellConn.*x/cellSize(orthDir),stabWeights,...
    'UniformOutput',false) ;
ccRowDiagW = cellfun(@(x) sparse(diag(sum(x,2))),cellConnW,...
    'UniformOutput',false) ;
ccColDiagW = cellfun(@(x) sparse(diag(sum(x,1))),cellConnW,...
    'UniformOutput',false) ;

%% Diffusion

Ks = getConductivity(assembler) ;
Korder2get = min(dir,order-1) ; % dir if order==3, 1 if order==2
Ks = mat2cell(Ks.space.spaces{Korder2get},Ks.space.sz(Korder2get),...
    ones(1,Ks.space.dim(Korder2get))) ; % cell arry is single row
IdK = cellfun(@(x) sparse(diag(x)),Ks,'UniformOutput',false) ;
% IdK==IdKW <= ~useAverageWeights

%% Consistency

if getUseAverageWeights(assembler)
    avgWeights = getAverageWeights(assembler,min(dir,order-1)) ;
    if ~iscell(avgWeights)
        avgWeights = {avgWeights} ;
    elseif size(avgWeights,1)~=1 % ensure single row format
        avgWeights = avgWeights(:)' ;
    end
else
    % Apply factor 1/2 (order 2) or 1/sqrt(2) (order 3)
    avgWeights = {ones(size(cellConn))/2^(1/(order-1))} ;
end

% For consistency operator, apply conductivity to connectivity matrix
cellConnK = cellfun(@(x) cellConn*x,IdK,'UniformOutput',false) ;
r_avgW = length(avgWeights) ;
r_K = length(IdK) ;
IdKW = cell(1,r_K*r_avgW) ;
cellConnKW = cell(1,r_K*r_avgW) ;
for n = 1:r_avgW
    ind = (1:r_K)+(n-1)*r_K ;
    IdKW(ind) = cellfun(@(x) x.*avgWeights{n},IdK,'UniformOutput',false) ;
    cellConnKW(ind) = cellfun(@(x) x.*avgWeights{n},cellConnK,...
        'UniformOutput',false) ;
end
ccRowDiagKW = cellfun(@(x) sparse(diag(sum(x,2))),cellConnKW,...
    'UniformOutput',false) ;
ccColDiagKW = cellfun(@(x) sparse(diag(sum(x,1))),cellConnKW,...
    'UniformOutput',false) ;

%% Mean penalization

if order == 3
    switch getConstantNullification(assembler) ;
        case '1point'
            cstNullOperator = {sparse(1,1,1,cellNb,cellNb)} ; % Lower Left cell
        case 'full'
            cstNullOperator = {eye(cellNb,cellNb)} ;
        otherwise
            error('QPDiffusionAssembler: unknown constNullification method')
    end
else % unused in order 2 structure (assembled in assembleMesoOperators)
    cstNullOperator = {} ;
end

%% Storage

% Operators are stored in columns of cell arrays. Each K_n^I or \omega_n^I
% is associated to a column. This results from weightOperator and Ks single
% row format.

diffusionOperators = IdK ;

consistencyOperators = [ccRowDiagKW ; ... % l({}^t\ki^q\odot K_n^I)
    ccColDiagKW ; ... % l(\ki^q\odot K_n^I)
    cellfun(@transpose,cellConnKW,'UniformOutput',false) ; ... % {}^t\ki^q\odot K_n^I
    cellConnKW] ; % \ki^q\odot K_n^I

stabilisationOperators = [ccRowDiagW ; ... % l({}^t\ki^q\odot{}^t\omega)
    ccColDiagW ; ... % l(\ki^q\odot\omega)
    cellfun(@transpose,cellConnW,'UniformOutput',false) ; ... % {}^t\ki^q\odot{}^t\omega
    cellConnW] ; % \ki^q\odot\omega

if order == 3 % else do not include (keep TSpace as small as possible)
    consistencyOperators = [consistencyOperators ; IdKW] ; % diag(K_n^I)
    stabilisationOperators = [stabilisationOperators ; IdW] ; % diag(\omega)
end

% Store as structure
mesoOperators = struct('diffusion',{diffusionOperators},'consistency',...
    {consistencyOperators},'stabilisation',{stabilisationOperators},...
    'cstNull',{cstNullOperator}) ;

end