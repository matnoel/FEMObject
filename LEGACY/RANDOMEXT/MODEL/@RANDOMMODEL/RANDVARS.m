function rvs = RANDVARS(S)
% function rvs = RANDVARS(S)

rvs = RANDVARS(RANDVARS(S.ls),RANDVARS(MATERIALS(S)));
