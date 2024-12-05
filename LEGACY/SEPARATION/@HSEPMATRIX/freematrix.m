function H = freematrix(H,var,model)
% function H = freematrix(H,var,model)
% H : HSEPMATRIX,
% var  : variable cible
% model : modele donnant l'info

if length(var)==1
    var2dim=getvar2dim(H.tree);
    path=var2dim{var};
else
    path=var;
end
D=path(1);
if isa(H.F{1,D},'HSEPMATRIX')
    path=path(2:end);
    for r=1:H.m
        H.F{r,D}=freematrix(H.F{r,D},path,model);
    end
elseif isa(H.F{1,D},'SEPMATRIX')
    for r=1:H.m
        H.F{r,D}=freematrix(H.F{r,D},path(2),model);
    end
end





