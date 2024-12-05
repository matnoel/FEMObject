function tsz = tensorSize(model)
% tsz = tensorSize(model)
% Provide tensor size in TSpaceVectors format.
% Intended for TuckerLikeTensor static constructors.

cellNum = getCellNum(model) ;
tsz = [formatIndex(getOrder(model),cellNum,cellNum), ...
    getNbCellDoF(model)] ;
end

% % Alternative
% switch getOrder(model)
%     case 2
%         tsz = [getCellNb(model) getNbCellDoF(model)] ;
%     case 3
%         tsz = [getCellNum(model) getNbCellDoF(model)] ;
% end