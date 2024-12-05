function dWeights = directionalWeights(weights,cellNum)

if isfloat(weights) ; weights = {weights} ; end

if numel(weights) == 1 % cell containing weight operator matrix
    sz = size(weights{1},1) ;
    if sz > 1 % safety
        dWeights = cell(4,1) ; % sort first by direction
        dWeights{1} = {[diag(weights{1},1) ; diag(weights{1},1-sz)]} ;
        dWeights{2} = {[diag(weights{1},cellNum(1)) ; ...
            diag(weights{1},cellNum(1)-sz)]} ;
        dWeights{3} = {[diag(weights{1},sz-1) ; diag(weights{1},-1)]} ;
        dWeights{4} = {[diag(weights{1},sz-cellNum(1)) ; ...
            diag(weights{1},-cellNum(1))]} ;
    else
        dWeights = repmat({weights},4,1) ;
    end
else
    dWeights = repmat(cell(2,1),4,1) ; % sort by direction then order
    sz1 = size(weights{1}{1},1) ;
    sz2 = size(weights{2}{1},1) ;
    for n = 1:numel(weights{1}) % numel(weights{1})==numel(weights{2})
        dWeights{1}{2} = cell2mat(cellfun(@diag, ...
            weights{2}(:)','UniformOutput',false)) ;
        dWeights{2}{1} = cell2mat(cellfun(@diag, ...
            weights{1}(:)','UniformOutput',false)) ;
        dWeights{3}{2} = dWeights{1}{2} ;
        dWeights{4}{1} = dWeights{2}{1} ;
        if sz1 > 1 % safety
            dWeights{1}{1} = cell2mat(cellfun(@(x) [diag(x,1) ; diag(x,1-sz1)],...
                weights{1}(:)','UniformOutput',false)) ;
            dWeights{3}{1} = cell2mat(cellfun(@(x) [diag(x,sz1-1) ; diag(x,-1)],...
                weights{1}(:)','UniformOutput',false)) ;
        else
            dWeights{1}{1} = dWeights{2}{1} ;
            dWeights{3}{1} = dWeights{2}{1} ;
        end
        if sz2 > 2 % safety
            dWeights{2}{2} = cell2mat(cellfun(@(x) [diag(x,1) ; diag(x,1-sz2)],...
                weights{2}(:)','UniformOutput',false)) ;
            dWeights{4}{2} = cell2mat(cellfun(@(x) [diag(x,sz2-1) ; diag(x,-1)],...
                weights{2}(:)','UniformOutput',false)) ;
        else
            dWeights{2}{2} = dWeights{1}{2} ;
            dWeights{4}{2} = dWeights{1}{2} ;
        end
    end
end
end