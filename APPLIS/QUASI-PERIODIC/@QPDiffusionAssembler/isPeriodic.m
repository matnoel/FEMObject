function bool = isPeriodic(assembler,direction)
% bool = isPeriodic(assembler,direction)
% Returns true if assembler has PBC along direction in second input 
% argument (1 and/or 2).

if nargin < 2
    direction = [1 2] ;
end

%% Recursion loop
if numel(direction) > 1
    bool = true ;
    for d = direction
        bool = bool && isPeriodic(assembler,d) ;
    end
    return ;
end

%% Method

bc = getBC(assembler) ;
bool = ischar(bc) && strcmp(bc,'PBC') ;

if iscell(bc)
    cellNum = getCellNum(assembler) ;
    pbcLoc = cellfun(@(x) getType(x)==5,bc) ;
    orthDir = setdiff([1 2],direction) ;
   for i = find(pbcLoc(:))'
       cells = getCells(bc{i}) ;
       cells = formatIndex(3,cellNum,cells{direction}) ;
       bool = all(cells(:,direction)==cellNum(direction))...
           && all(sort(cells(:,orthDir)) ==  (1:cellNum(orthDir))') ;
       if bool ; return ; end
   end
end

end