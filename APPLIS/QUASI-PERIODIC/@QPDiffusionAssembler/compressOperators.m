function [assembler,compressionTime] = compressOperators(assembler,tolerance)
% [assembler, compressionTime] = compressOperators(assembler,tolerance)

if nargin == 1
    tolerance = getTolSVD(assembler) ;
end

compressionClock = tic ;

compressor = Truncator('tolerance',tolerance) ;

% Compress SWIP and source operators
lhsOp = compressor.truncate(getLHSOperator(assembler)) ;
rhsOp = compressor.truncate(getRHSOperator(assembler)) ;

% Convert and store
[lhsOp,rhsOp] = convertTensors(lhsOp,rhsOp) ;
assembler = setLHSOperator(assembler,lhsOp) ;
assembler = setRHSOperator(assembler,rhsOp) ;

compressionTime = toc(compressionClock) ;
end