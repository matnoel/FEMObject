function assembler = assemblePatterns(assembler,patternsTable,patterns)
% Assembler = assemblePatterns(assembler,patternsTable,patterns)

if nargin < 3
    patterns = getPatterns(assembler) ;
    if nargin < 2
        patternsTable = getPatternsTable(assembler) ;
    end
end

patternFields = drawCellPattern(getCellCoord(assembler),patterns) ;
cellFields = patternFields*patternsTable ;
assembler = setFields(assembler,cellFields) ;

end