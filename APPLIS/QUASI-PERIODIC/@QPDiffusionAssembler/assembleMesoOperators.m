function [Assembler,assemblerTime] = assembleMesoOperators(Assembler)
% [Assembler,assemblerTime] = assembleMesoOperators(Assembler)

assemblerClock = tic ;

%% Diffusion, consistency, stabilisation (and cstNull if order 3) operators
mesoOperators = assembleMesoOperatorsDirection(Assembler,1) ;
newOps = assembleMesoOperatorsDirection(Assembler,2) ;
opNames = fieldnames(mesoOperators) ;
for i = 1:numel(opNames)
    opName = opNames{i} ;
    switch getOrder(Assembler)
        case 2
            % Avoid non-directional diffusion operators redundancy if order = 2
            if strcmp(opName,'diffusion')
                mesoOperators.(opName) = {mesoOperators.(opName)} ;
            else
                mesoOperators.(opName) = {[mesoOperators.(opName) ; newOps.(opName)]} ;
            end
        case 3
            mesoOperators.(opName) = {mesoOperators.(opName) ; newOps.(opName)} ;
    end
end

%% Mean penalization

if getOrder(Assembler) == 2 % else it has been done in assembleOrderMesoOperators
    cellNb = getCellNb(Assembler) ;
    switch getConstantNullification(Assembler) ;
        case '1point'
            mesoOperators.cstNull = {{sparse(1,1,1,cellNb,cellNb)}} ; % Lower Left cell
        case 'full'
            mesoOperators.cstNull = {{eye(cellNb,cellNb)}} ;
        otherwise
            error('QPDiffusionAssembler: unknown constNullification method')
    end
end

%% Source

source = getSource(Assembler) ;
if ~ischar(source) && isa(source.space,'TSpaceVectors')
    mesoOperators.rhs = source.space.spaces(1:end-1) ;
else
    mesoOperators.rhs = [] ;
end

%% Storage

Assembler = setMesoOperators(Assembler,mesoOperators) ;

assemblerTime = toc(assemblerClock) ;
end