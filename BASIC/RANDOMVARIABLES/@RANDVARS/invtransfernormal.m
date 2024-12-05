function ax = invtransfernormal(a,xix)
% function ax = invtransfernormal(a,xix)
% a : RANDVARS contenant M variables aleatoires   
% variables supposees independantes
% xix, ax : vecteur de taille n*M
% ax(:,i)=F_ai^-1(Phi(xix(:,i)))
ax=zeros(size(xix));
a=getrandvar(a);
for i=1:size(ax,2)
    ax(:,i)=invtransfernormal(a{i},xix(:,i));
end