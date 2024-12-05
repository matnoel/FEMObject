function distribution = distribute(assembler)
% distribution = distribute(assembler)

distributorFunc = getDistributor(assembler) ;
cellNb = getCellNb(assembler) ;
distribution = distributorFunc(cellNb) ;

end