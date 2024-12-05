function rep = polycmp(PC1,PC2)
% function rep = polycmp(PC1,PC2)

rep = polycmp(PC1.RANDPOLYS,PC2.RANDPOLYS);
if rep
    rep = rep & (PC1.typebase == PC2.typebase);
    rep = rep & (~any(PC1.p ~= PC2.p));
    rep = rep & (~any(PC1.n ~= PC2.n));
    rep = rep & (~any(PC1.P ~= PC2.P));
end

