function indicator = mesoIndicator(model,indices,isOperator)
% indicator = mesoIndicator(model,indices,isOperator)
% Convert cells indices into characteristics operators (or indicators)

if nargin<3
    isOperator = false ;
end

assert(isnumeric(indices)||isnumeric(indices{1}),...
    'Second input must be double or cell array of double')

%% Cell recursion

if iscell(indices)
    cNb = numel(indices) ;
    order = getOrder(model) ;
    indicator = cell(order-1,cNb) ;
    for c = 1:cNb
        indicator(:,c) = mesoIndicator(model,indices{c},isOperator) ;
    end
    return
end

% Safety
if isempty(indices)
    indicator = cell(1,getOrder(model)-1) ;
    return
end

%% Conversion from indices to indicators/operators

% Pre-processing
cellNum = getCellNum(model) ;
cellNb = prod(cellNum) ;
order = getOrder(model) ;

% Build indicators/operators
if isOperator
    if size(indices,2)>2
        indices(:,1) = formatIndex(2,cellNum,indices(:,1:end/2)) ;
        indices(:,2) = formatIndex(2,cellNum,indices(:,1+end/2:end)) ;
    end
    % Matrix with 1 on elements in indices, 0 everywhere else
    rep = full(sparse(indices(:,1),indices(:,2),1,cellNb,cellNb)) ;
else % then is vector
    indices = formatIndex(2,cellNum,indices) ;
    % Column with 1 on lines in indices, 0 everywhere else.
    rep = full(sparse(indices,1,1,cellNb,1)) ;
end

if order == 2 % Stop here for order 2.
    indicator = {rep} ; % Store in TSpaceVectors format
    if isOperator
        indicator = {indicator} ; % Store in TSpaceOperators format
    end
else
    indicator = svdQP(model,rep) ; % Order 3 processing (truncated SVD)
    % Already formatted
end


end