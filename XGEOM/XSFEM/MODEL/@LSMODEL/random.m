function [S,r] = random(S)
% function [S,r] = random(S)

RV = RANDVARS(S);
r = random(RV);
S = randomeval(S,r,RV);
% S = lssplitelem(S);
