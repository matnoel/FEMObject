function t = mesoIndicatorT(model,indices)
% indicator = mesoIndicatorT(model,indices)
% Convert cells indices into characteristics indicators.
% Indices are assumed to be correctly formatted, as with
% formatIndex(order,cellNum,indices).
% Input indices with double number of column will return operators.

% Safety
if isempty(indices)
    t = TuckerLikeTensor.zeros(tensorSize(model)) ;
    return
end

%% Conversion from indices to indicators

% Pre-processing
cellNum = getCellNum(model) ;
cellNb = prod(cellNum) ;
order = getOrder(model) ;
isOperator = size(indices,2) == 2*(order-1) ;

% Build indicators for order 2 format
if isOperator % operators
    if order == 3
        indices = [formatIndex(2,cellNum,indices(:,1:end/2)) ...
            formatIndex(2,cellNum,indices(:,1+end/2:end))] ;
    end
    % Matrix with 1 on elements in indices, 0 everywhere else
    rep = sparse(indices(:,1),indices(:,2),1,cellNb,cellNb) ;
else % vectors
    indices = formatIndex(2,cellNum,indices) ;
    % Column with 1 on lines in indices, 0 everywhere else.
    rep = sparse(indices,1,1,cellNb,1) ;
end

if order == 2 % Stop here for order 2.
    indicator = {rep} ; % Store in TSpaceVectors format
    if isOperator
        indicator = {indicator} ; % Store in TSpaceOperators format
    end
else % convert for order 3
    indicator = svdQP(model,rep) ; % Order 3 processing (truncated SVD)
    % Already formatted, unless degenerate case (cf. infra)
    if isOperator && isscalar(rep)
       indicator = cellfun(@(x){x},indicator,'UniformOutput',false) ;
    end
end

%% Conversion from indicator to tensor

% TSpace
space = cell(order,1) ;
space(1:order-1) = indicator ;
if isOperator
%     space{end} = {ones(getNbCellDoF(model))} ;
    space{end} = {speye(getNbCellDoF(model))} ;
    space = TSpaceOperators(space) ;
else
    space{end} = ones(getNbCellDoF(model),1) ;
    space = TSpaceVectors(space) ;
end

% Core
if all(space.dim==1)
    core = DiagonalTensor(ones(min(space.dim),1),order) ;
else % order 3 where space.dim is [d d 1] with d>1
    core = zeros(space.dim) ;
    core(:,:,1) = eye(space.dim(1:end-1)) ;
    core = FullTensor(core,order,space.dim) ;
end

t = TuckerLikeTensor(core,space) ;
end