function S = LSMODEL(M,ls)

if nargin==0
    S = struct();
    S.ls = LEVELSETS();
    S.H = [];
    S = class(S,'LSMODEL',MODEL());
    S = createlsmodelfaces(S);

elseif nargin==1 && isa(M,'LSMODEL')
    S = M;
elseif nargin==1 && isa(M,'MODEL')
    S = LSMODEL(M,LEVELSETS());
    S = createlsmodelfaces(S);

else

    S = struct();
    if nargin==1
        S.ls = LEVELSETS();
    else
        S.ls = LEVELSETS(ls);
    end
    S.H = [];
    S = class(S,'LSMODEL',getmodel(M));
    S = createlsmodelfaces(S);

end
