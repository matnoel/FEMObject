function [pb,time] = updateSource(pb,source)
% [pb,time] = updateSource(pb,source)

if nargin < 2
    source = getSource(pb) ;
end

[assembler,time] = updateSource(getOperatorAssembler(pb),source) ;
pb = updateOperatorAssembler(pb,assembler) ;

end