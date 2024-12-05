function [assembler,time] = updateBC(assembler,bc,bcNum)
% [assembler,time] = updateBC(assembler,bc,bcNum)

if nargin < 3
    bcNum = [] ;
    if nargin < 2
        bc = getBC(assembler) ;
    end
end

clock = tic ;

if ~iscell(bc)
    bc = {bc} ;
end

newBC = getBC(assembler) ;
if ~isempty(bcNum)
    newBC(bcNum) = bc ;
else
    newBC = bc ;
end

assembler = setBC(assembler,newBC) ;
assembler = assemble(assembler,[]) ;

time = toc(clock) ;
end