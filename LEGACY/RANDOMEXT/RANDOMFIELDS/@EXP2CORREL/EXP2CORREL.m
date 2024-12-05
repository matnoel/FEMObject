function co = EXP2CORREL(L)
% function co = EXP2CORREL(L)

co.L=L;
cop = RFCORREL('exp2',co);
co = class(co,'EXP2CORREL',cop);
