function [assembler,time] = updateSource(assembler,source)
% [assembler,time] = updateSource(assembler,source)

if nargin < 2
    source = getSource(assembler) ;
else
    source = formatSource(assembler,source) ;
end

assembler = setSource(assembler,source) ;
[assembler,time] = assemble(assembler) ;

end