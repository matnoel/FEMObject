function [assembler,assemblerTime] = assemble(assembler)
% [Assembler,assemblerTime] = assemble(Assembler)

assemblerClock = tic ;

% For simplicity and speed
order = getOrder(assembler) ;
dist = getDistribution(assembler) ;
fields = getFields(assembler) ;
model = getModel(assembler) ;

if isempty(fields)
    assembler = assemblePatterns(assembler) ;
    fields = getFields(assembler) ;
end

% Remove absent phases
phasesPresent = ~cellfun(@isempty,dist) ;
dist = dist(phasesPresent) ;
fields = fields(:,phasesPresent) ;

% Build space
Kspace = cell(order,1) ;
Kspace{end} = fields ;
mesoInd = mesoIndicator(model,dist,0) ;
for o = 1:order-1
    Kspace{o} = [mesoInd{o,:}] ; % store indicators as columns in matrix
end
Kspace = TSpaceVectors(Kspace) ;

% Build core
if order == 2
    Kcore = DiagonalTensor(ones(Kspace.dim(1),1),order) ;% each indicator is elementary
else % order 3
    coreSize = Kspace.dim ;
    indicatorRank = cellfun(@(x) size(x,2),mesoInd(1,:)) ;
    Kcore = zeros(coreSize) ; 
    for i = 1:coreSize(end)
        currentRank = 1+sum(indicatorRank(1:i-1)) ;
        range = currentRank:(currentRank+indicatorRank(i)-1) ;
        Kcore(range,range,i) = eye(indicatorRank(i)) ;
    end
    Kcore = FullTensor(Kcore,order,coreSize) ;
end

K = TuckerLikeTensor(Kcore,Kspace) ;
assembler = setConductivity(assembler,K) ;
assembler = setConductivityBounds(assembler) ;

assemblerTime = toc(assemblerClock) ;

end