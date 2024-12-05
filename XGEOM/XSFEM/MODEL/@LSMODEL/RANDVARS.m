function rvs = RANDVARS(S)
% function rvs = RANDVARS(S)

rvs = RANDVARS(RANDVARS(S.MODEL),RANDVARS(S.ls));
