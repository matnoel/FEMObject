function H = unfreevector(H,var,model)
% function H = unfreevector(H,var,model)
% H : HSEPMATRIX,
% var  : variable cible
% model : modele donnant l'info

if length(var)==1
    T=H.tree;
    var2dim=getCvar(T);
    path=var2dim{var};
end
D=path(1);


if isa(H.F{1,D},'HSEPMATRIX')
    path=path(2:end);
    for r=1:H.m
        H.F{r,D}=unfreevector(H.F{r,D},path,model);
    end
elseif isa(H.F{1,D},'SEPMATRIX')
    for r=1:H.m
        H.F{r,D}=unfreevector(H.F{r,D},path(2),model);
    end
end


