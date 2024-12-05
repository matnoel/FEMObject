function lam=thermique_nonlin_GM(U,uM,tol,S,PC,X)

if ~isempty(uM)
    M = getm(uM);
else
    M=0;    
end

thermique_nonlin_forms

if M==0
    a1 = X{1}*a{S}(U,U);
    a2=0;
    a3 = X{2}*n{S}(U,U);
    if nbsec==1
        b = X{3}*l{S}(U);
    else
        b = X{3}*l{S}(U)+X{4}*l2{S}(U);        
    end
else
    stoU=double(getV(uM));
    stol=getL(uM);


    N1 = g{S}(U,U,:,:)+2*g{S}(U,:,U,:);
    a1 = X{1}*a{S}(U,U) + X{2}*stol'*(stoU'*N1*stoU)*stol;
    N2 = g{S}(U,:,U,U)+2*g{S}(U,U,:,U);
    a2 = (X{2}*stol)'*(stoU'*N2);
    a3 = X{2}*n{S}(U,U);
    if nbsec==1
        b = X{3}*l{S}(U);
    else
        b = X{3}*l{S}(U)+X{4}*l2{S}(U);    
    end

    b = b- X{1}*(stol'*(stoU'*a{S}(:,U)));
    for j=1:getm(uM)
        b = b - (X{2}*stol(j))*(stol'*(stoU'*g{S}(U,stoU(:,j),:,:)*stoU)*stol);
    end
end

lam = GSDthermiquesolvesto(PC,tol,a1,a2,a3,b);



