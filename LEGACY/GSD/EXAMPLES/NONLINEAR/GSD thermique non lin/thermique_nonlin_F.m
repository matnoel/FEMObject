function W=thermique_nonlin_F(Lam,tol,S,X)
thermique_nonlin_forms
s2=numel(Lam);
ka=expectmtimes(X{1},Lam,Lam');
kn = zeros(s2,s2,s2,s2);
for k=1:s2
    for j=1:s2
        kn(:,:,k,j)=expectmtimes(X{2}*Lam,Lam',Lam(k)*Lam(j));
        kn(:,:,j,k)=kn(:,:,k,j);
    end
end
if nbsec==1
    kl=expectmtimes(X{3},Lam);
else
    kl{1}=expectmtimes(X{3},Lam);
    kl{2}=expectmtimes(X{4},Lam);
end

W = GSDthermiquesolvedetreac(S,tol,ka,kn,kl);




