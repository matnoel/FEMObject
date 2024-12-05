function [assembler,assemblerTime] = assemble(assembler,orders2Assemble)
% [Assembler,assemblerTime] = assemble(Assembler,orders2Assemble)

order = getOrder(assembler);

if nargin == 1
    orders2Assemble = 1:order ;
end

KTime = 0 ;
mesoTime = 0 ;
microTime = 0 ;
compTime = 0 ;
patchTime = 0 ;

if isempty(getConductivity(assembler))
    cA = getConductivityAssembler(assembler) ;
    [cA,KTime] = assemble(cA) ;
    tic ;
    assembler = updateConductivityAssembler(assembler,cA) ;
    assemblerTime = KTime + toc ;
    return % update method calls assemble method
end

if any(ismember(1:order-1,orders2Assemble)) || KTime > 0
    if ~all(ismember(1:order-1,orders2Assemble))
        warning('Forcing complete meso assembling anyway')
    end
    [assembler,mesoTime] = assembleMesoOperators(assembler) ;
end

if ismember(order,orders2Assemble) || KTime > 0
    if ~ismember(order,orders2Assemble)
        warning('Forcing micro assembling anyway')
    end
    [assembler,microTime] = assembleMicroOperators(assembler) ;
end

% Assembling operators
[diffOp,diffTime] = assembleDiffusionOperator(assembler) ;
[consOp,consTime] = assembleConsistencyOperator(assembler) ;
[stabOp,stabTime] = assembleStabilisationOperator(assembler) ;
% [cstNullOp,cstNullTime] = assembleCstNullOperator(assembler) ;
[sourceOp,sourceTime] = assembleSourceOperator(assembler) ;
[lhsBCOp,rhsBCOp,bcTime] = assembleBCOperator(assembler) ;
ps = getPatches(assembler) ;
if ~isempty(ps)
    [ps,patchTime] = assemble(ps,getMicroOperators(assembler),...
        getPenalty(assembler)) ;
    assembler = setPatches(assembler,ps) ;
end

% Add, convert and store
lhsOp = diffOp - consOp - consOp' + stabOp + lhsBCOp ;
rhsOp = sourceOp + rhsBCOp ;
[lhsOp,rhsOp] = convertTensors(lhsOp,rhsOp) ;
assembler = setLHSOperator(assembler,lhsOp) ;
assembler = setRHSOperator(assembler,rhsOp) ;

% Compression
if getUseCompression(assembler)
    [assembler,compTime] = compressOperators(assembler) ;
end

assemblerTime = KTime + mesoTime + microTime + diffTime + consTime ...
    + stabTime + sourceTime + bcTime + compTime + patchTime ;
end