function ps = updateConductivityAssembler(ps,cA)
% ps = updateConductivityAssembler(ps,cA)

if nargin == 1
    cA = getConductivityAssembler(ps) ;
end

cA = realConductivityAssembler(ps,cA,1) ;
ps = setConductivityAssembler(ps,cA) ;
model = getModel(ps) ;
pCells = getCells(ps) ;
patches = getPatches(ps) ;
for i = 1:size(pCells,1)
    pModel = subModel(model,pCells{i}{1}) ;
    patches{i} = setModel(patches{i},pModel) ;
end

end