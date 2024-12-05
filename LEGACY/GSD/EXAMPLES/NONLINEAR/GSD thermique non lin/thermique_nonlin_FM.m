function U=thermique_nonlin_FM(lam,uM,tol,S,X)
thermique_nonlin_forms

ka=expectmtimes(X{1},lam,lam);
lam2 = lam*lam;
kn=expectmtimes(X{2},lam2,lam2);
if nbsec==1
    kl=expectmtimes(X{3},lam);
else
    kl{1}=expectmtimes(X{3},lam);
    kl{2}=expectmtimes(X{4},lam);
end

if ~isempty(uM)
    M = getm(uM);
else
    M=0;    
end

if M==0
    U = GSDthermiquesolvedet(S,tol,ka,kn,kl);
else
    stoU=double(getV(uM));
    stol=getL(uM);
    stol = calc_masse(stol);
    lamstol = lam*stol;
    Utilde = stoU * expectmtimes(lamstol,X{2},lam2);
    Urond =  stoU * expectmtimes(stol,X{1},lam);
    Ubar = cell(1,getm(uM));
    Uchap = cell(getm(uM),getm(uM));
    Ui = cell(1,getm(uM));
    for j=1:getm(uM)
        Ui{j}= stoU(:,j);
        Ubar{j}=stoU * expectmtimes(stol,stol(j),X{2}*lam2);
        for k=j:getm(uM)
            Uchap{j,k} = stoU * expectmtimes(X{2},stol(j)*stol(k),lamstol);
            Uchap{k,j} = Uchap{j,k};
        end
    end
    U = GSDthermiquesolvedetM(S,tol,ka,kn,kl,Ui,Urond,Utilde,Ubar,Uchap);

end


