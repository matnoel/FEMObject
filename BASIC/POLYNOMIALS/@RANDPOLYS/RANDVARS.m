function rvs = RANDVARS(H)
% function rvs = RANDVARS(H)

H = RANDPOLYS(H);
rv = cell(1,H.M);
for k=1:H.M
    rv{k}=RANDVAR(H.h{k});
end

rvs = RANDVARS(rv{:});
