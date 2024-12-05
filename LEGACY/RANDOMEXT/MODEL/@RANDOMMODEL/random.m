function S = random(S)
% function S = random(S)

RV = RANDVARS(S);
r = random(RV);
S = randomeval(S,r,RV);
S = lssplitelem(S);
