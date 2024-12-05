function co = EXPCORREL(L)
% function co = EXPCORREL(L)

co.L=L;
cop = RFCORREL('exp',co);
co = class(co,'EXPCORREL',cop);
