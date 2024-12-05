function gauss=calc_smolyak(H,l,varargin)

% function gauss=calc_gausspoints(H,l)
%
% calcul d'une quadrature de Smolyak
% l : nombre de points de Gauss maxi dans les quadratures
% unidimensionnelles (niveau maxi de la quadrature de smolyak)
% si l est un vecteur de taille z : l(j) donne le nombre de points de Gauss
% de la quadrature unidimensionnelles d'indice i

m=H.M;

if length(l)>1
    Q = cell(m,length(l));
for i=1:m
    for j=1:l
Q{i,j}=calc_gausspoints(H.h{i},l(j));
    end
end
l=length(l);
else
    Q = cell(m,l);
for i=1:m
    for j=1:l
Q{i,j}=calc_gausspoints(H.h{i},j);
    end
end
end

PC = POLYCHAOS(m,l+m-1);
ind = getindices(PC);
ind = ind(:,1:end-1);
rep =  find(all(ind>0,2) & sum(ind,2)>=l & sum(ind,2)<=l+m-1);
ind = ind(rep,:);

gauss.coord = zeros(0,m);
gauss.w = zeros(1,0);
gauss1D=cell(1,m);

if ischarin('count',varargin)
gauss=0;
for k=1:size(ind,1)
    
    nn=1;
    for i=1:m
        nn=nn*Q{i,ind(k,i)}.nbgauss;
    end
    gauss = gauss+nn;
end
return
    
end

for k=1:size(ind,1)
for i=1:m
    gauss1D{i} = Q{i,ind(k,i)};
end
    gaussmD = tensorize_quadrature_rule(gauss1D);
    gaussmD.w = gaussmD.w*(-1)^(m+l-1-sum(ind(k,:)))*nchoosek(m-1,sum(ind(k,:))-l);  
%fprintf('---\n indices %d , nombre de points = %d-----',k,length(gaussmD.w))
%disp(ind(k,:))
    gauss.coord = [gauss.coord;gaussmD.coord];
    gauss.w = [gauss.w,gaussmD.w];    
end
gauss.nbgauss = length(gauss.w);
%gauss.coord

if ischarin('unique',varargin)
prec = 1e-8;
[B,I,J] = unique(floor(gauss.coord/prec)*prec,'rows');
fprintf('Nombre de points avant elimination  = %d\n',gauss.nbgauss);

if size(B,1)<gauss.nbgauss
gauss.coord=B;
w = gauss.w;
gauss.w = gauss.w(I);
K = setdiff(1:gauss.nbgauss,I);
for i=1:length(K)
gauss.w(J(K(i)))=gauss.w(J(K(i)))+w(K(i));    
end

gauss.nbgauss = length(gauss.w);
end
fprintf('Nombre de points apres elimination  = %d\n',gauss.nbgauss);
end

