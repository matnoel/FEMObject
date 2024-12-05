function f = untensorize(model,t,toSmooth)
% f = untensorize(model,t,toSmooth)

if nargin < 3
    toSmooth = true ;
end

%% Cell Recursion

if iscell(t)
    f = cell(size(t)) ;
    for i = 1:numel(t)
        f{i} = untensorize(model,t{i}) ;
    end
    return
end

%% Method

if isempty(t) ; f = [] ; return ; end % safety

f = doubleQP(t) ;

if toSmooth
    f = smooth(model,f) ;
end
end

%% Previous method (bugged for operators)

% % Method
% 
% tOrder = getOrder(model) ;
% isOperator = isa(t.space,'TSpaceOperators') ;
% 
% if isOperator
%     t = vectorize(t) ;
% end
% 
% 
% mesoDim = t.space.sz(1,1:tOrder-1) ;
% mesoDoFNb = prod(mesoDim) ;
% nbCellDoF = t.space.sz(end) ;
% 
% cellsInd = formatIndex(tOrder,mesoDim,(1:mesoDoFNb)') ;
% 
% cellValues = evalAtIndices(t,cellsInd,1:tOrder-1) ;
% cellValues = double(cellValues)' ; % one column per cell
% 
% if isOperator
%     f = zeros(sqrt(mesoDoFNb*nbCellDoF)) ;
%     cellsInd = formatIndex(2,mesoDim,cellsInd) ; % ensure indices
%     for c = 1:mesoDoFNb
%         cellOp = reshape(cellValues(:,c),sqrt(nbCellDoF),sqrt(nbCellDoF)) ;
%         cellPlace = zeros(sqrt(mesoDoFNb)) ;
%         cellPlace(cellsInd(c)) = 1 ;
%         f = f + kron(cellPlace,cellOp) ;
%     end
% else
%     f = cellValues(:) ;
% end
% 
% if toSmooth
%     f = smooth(model,f) ;
% end