function S=HSM2SM(H)
% S est l'expand de H en SEPMATRIX

totrank=totalrank(H);
totdim =getdim(H.tree);
F=cell(totrank,totdim);
alpha = ones(totrank,H.dim+1);
R=1;
for r=1:H.m
    % Quelques preparatifs : le rang de chaque branches
    TRD=zeros(H.dim,1);
    for d=1:H.dim
        if isa(H.F{r,d},'HSEPMATRIX')
            TRD(d) = totalrank(H.F{r,d});
        else
            TRD(d) = getm(H.F{r,d});
        end
    end
    trr=prod(TRD);
    m=1;
    M=trr;
    
    % On distribue les vecteurs aux bons endroits
    D=1;
    for d=1:H.dim
        if isa(H.F{r,d},'HSEPMATRIX')
            HS = HSM2SM(H.F{r,d});
        else
            HS = H.F{r,d};
        end
        rr=HS.m;
        dd=HS.dim;
        M=M/rr;
        F([R : R-1+trr],[D : D-1+dd]) = repmat(reshape(repmat(reshape(HS.F,1,rr*dd),m,1) ,rr*m,dd),M,1);
        alpha([R : R-1+trr],d) = repmat(reshape(repmat(reshape(HS.alpha,1,rr),m,1) ,rr*m,1),M,1);
        D=D+dd;
        m=m*rr;
    end
    alpha([R : R-1+trr],H.dim+1)=ones(trr,1)*H.alpha(r);
    R=R+trr;
end

alpha=prod(alpha,2)';

S=SEPMATRIX(F);
S.alpha=alpha;

% S=normalizefuns(S);
% [S.alpha pass]=sort(S.alpha,'descend');
% S.F(:,:)=S.F(pass,:);




