function rvpc = PCRANDVAR(u)
if isa(u,'PCRANDVAR')
    rvpc=u;
else

    if ~isa(u,'PCMATRIX') || numel(u)>1
        error('input must be a PCMATRIX of length 1')
    end

    [domain , domainborne] = domain(u);
    param = struct();
    rvpc = struct();
    rvp = RANDVAR('pc',param,domain,domainborne);
    rvpc = class(rvpc,'PCRANDVAR',u,rvp);
end

