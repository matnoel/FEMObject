function [flag,loc] = checkCells(ps)
% flag = checkCells(ps)
% Check lists of cells for every patch. Stops at first error and return
% corresponding flag.
% _loc = [patchNum listNum] may provide patch and (in some cases) list
% numbers of first error.
% _flag values:
%   1: correct.
%   0: overlap detected.
%   -1: cell number mismatch.
%   -2: nonconvex cell list detected.

loc = [] ;

% Check overlap
flag = checkOverlap(ps) ;
if ~flag ; return ; end

% Check number of cells and convexity
pCells = getCells(ps) ;
cellNum = getCellNum(ps) ;
pCells = formatIndex(2,cellNum,pCells) ;
for p=1:numel(pCells)
    % Check numbers of cells for current patch
    cellsNb = cellfun(@numel,pCells{p}) ;
    if any(diff(cellsNb)) % Cells numbers mismatch
        flag = -1 ;
        loc = p ;
        return
    % Check convexity
    elseif isprime(cellsNb(1)) && cellsNb(1)>2 % Cannot be convex : stop here
        flag = -2 ;
        loc = p ;
        return
    end
    for c=1:numel(pCells{p}) % Thorough test
        % We know there are no overlap so it is enough to check if cells
        % number is as expected for a rectangle spanning from min to max.
        subs = formatIndex(3,cellNum,pCells{p}{c}) ;
        expCellsNb = prod([1 1]+max(subs)-min(subs)) ; % expected value
        if expCellsNb ~= cellsNb(1)
            flag = -2 ;
            loc = [p c] ;
            return
        end
    end
end
end

