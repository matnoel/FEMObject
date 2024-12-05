function x = restrictTensor(model,x,cells,orthogonalise)
% x = restrictTensor(model,x,cells,orthogonalise)
% [Overload] calls tool function restrictTensor.

if nargin < 4
    orthogonalise = false ;
end

x = restrictTensor(x,cells,getCellNum(model),orthogonalise) ;

end

%% Previous method with subTensor (heavier)

% if isa(tensor.space,'TSpaceOperators')
%     
%     % Adapt meso coordinates for vectorized operators
%     csz2 = size(cells,2) ;
%     cellNum = getCellNum(model) ;
%     cellNb = getCellNb(model) ;
%     % Ensure two columns of cells indices
%     cells = [formatIndex(2,cellNum,cells(:,1:csz2/2)), ...
%         formatIndex(2,cellNum,cells(:,(1+csz2/2):csz2))] ;
%     % Consider as operator subscripts and apply format
%     cells = formatIndex(getOrder(model),cellNb*[1 1],cells) ;
%     
%     % Recurse on subTensor
%     tsz = tensor.sz ;
%     tensor = vectorize(tensor) ;
%     sTensor = subTensor(model,tensor,cells) ;
%     sTensor = unvectorize(sTensor,tsz) ;
%     return
% end
% 
% order = getOrder(model) ;
% cells = formatIndex(order,getCellNum(model),cells) ;
% 
% opValues = evalAtIndices(tensor,cells,1:order-1) ;
% opValues = double(opValues)' ;
% space = cell(order,1) ;
% space{order} = mat2cell(opValues,size(opValues,1),...
%     ones(1,size(opValues,2))) ;
% indic = mesoIndicator(model,cells,0) ;
% for o = 1:order-1
%     space{o} = [indic{o,:}] ; % store indicators as columns in matrix
% end
% space = TSpaceVectors(space) ;
% core = DiagonalTensor(ones(numel(opValues),1),order) ;
% sTensor = TuckerLikeTensor(core,space) ;