function w=SEPtimesSEP_FH_MULTM(u,v,dim)
if nargin<3, dim=1; end
if isa(u,'function_handle')
    w=u(v);
elseif isa(u,'MULTIMATRIX')
    w=getmatrix(switchmulti(v'*switchmulti(u)));
elseif isa(u,'SPARSETENSOR')
    w=SPTENSORtimesVECTOR(u,v,dim);
end

