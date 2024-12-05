function b=thermique_nonlin_res(uM,S,X,tol)

if nargin==0 || isempty(uM)
    M=0;
else
    if isa(uM,'PCMATRIX')
        if nargin<4
            tol=1e-8;
        end
        uM = spectral_decomposition(uM,'tol',tol);     
    end

    M = getm(uM);
    stoU=double(getV(uM));
    stol=getL(uM);

end


thermique_nonlin_forms

if nbsec==1
    b = X{3}*l{S}(:);
else
    b = X{3}*l{S}(:)+ X{4}*l2{S}(:);    
end

if M>0
    b = b - (a{S}(:,:)*stoU)*(X{1}*stol);
    for k=1:M
        for j=1:M
            temp = PCRADIALMATRIX(g{S}(:,:,stoU(:,k),stoU(:,j))*stoU,size(uM),stol);        
            b = b - temp*(X{2}*stol(k)*stol(j));
        end
    end
end
