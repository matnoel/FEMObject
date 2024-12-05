function bases = FunctionalBases(PC)
% function bases = FunctionalBases(PC)
% convert an object of type POLYCHAOSTP into a FunctionalBases
% PCTP : object of type POLYCHAOSTP
% bases : FunctionalBases

bases = cellfun(@PolychaosFunctionalBasis,getpcgroups(PC),...
                          'UniformOutput',false);
bases=FunctionalBases(bases);                 