function xix = transfernormal(a,ax)
% function xix = transfernormal(a,ax)
% a : RANDVARS contenant M variables aleatoires   
% variables supposees independantes
% xix, ax : vecteur de taille n*M
% xix(:,i)=Phi^-1(F_ai(ax(:,i)))
xix=zeros(size(ax));
a=getrandvar(a);
for i=1:size(ax,2)
    xix(:,i)=transfernormal(a{i},ax(:,i));
end