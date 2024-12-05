function [pb,time] = updateBC(pb,bc,bcNum)
% [pb,time] = updateBC(pb,bc,bcNum)

if nargin < 3
    bcNum = [] ;
    if nargin < 2
        bc = getBC(pb) ;
    end
end

[assembler,time] = updateBC(getOperatorAssembler(pb),bc,bcNum) ;
pb = updateOperatorAssembler(pb,assembler) ;

end