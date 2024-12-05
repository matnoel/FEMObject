function b = calc_nonhomogeneous_vector(BC,K)
% function b = calc_nonhomogeneous_vector(BC,K)

if all(gettypes(BC)==0) && ~all(ishomogeneous(BC))
    rep1 = BC.ddlfree;
    rep2 = BC.ddlbloque;
    v = zeros(length(rep1),1);
    v = unfreevector(BC,v);
    b = K(rep1,rep2)*v(rep2);
elseif ~all(ishomogeneous(BC))
    error('pas programme')
else
    b = zeros(getnbddlfree(BC),1);
end