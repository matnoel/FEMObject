function cells = formatCells(ps,cells)
% cells = formatCells(ps,cells)
% Format all cells indices cells for ps. Calls formatIndex.

% Safety
if isempty(cells) || (iscell(cells) && all(cellfun(@isempty,cells)))
    return
end

order = getOrder(ps) ;
cellNum = getCellNum(ps) ;
if iscell(cells{1})
    for p = 1:numel(cells)
        for c = 1:numel(cells{p})
            cells{p}{c} = formatIndex(order,cellNum,cells{p}{c}) ;
        end
    end
else
    for c = 1:numel(cells)
        cells{c} = formatIndex(order,cellNum,cells{c}) ;
    end
end
end

