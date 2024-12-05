function [gauss,Q]=calc_smolyak_curtis(H,l,base,varargin)

% function gauss=calc_smolyak_curtis(H,l,base,varargin)
%
% calcul d'une quadrature de Smolyak
% base : rapport du nombre de points entre niveaux
% l : niveau maxi de la quadrature (nombre de points base^(l-1)+1 dans la 
% quadrature fine)

m=H.M;
if nargin<=2 || isempty(base)  
    base=2;
end
if base<=1
    error('le nombre de points dans la petite grille doit etre superieur a 2')
end

   ni(1) = 1;
   for j=2:l
    ni(j) = 1+base^(j-1);
   end 
   
l=length(ni);


Q = cell(m,l);
for i=1:m
    for j=1:l
[Q{i,j}.coord,Q{i,j}.w]=fclencurt(ni(j),-1,1);
Q{i,j}.w=Q{i,j}.w';
Q{i,j}.nbgauss = ni(j);
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
    prec = 1e-12;
[B,I,J] = unique(floor(gauss.coord/prec)*prec,'rows');

if size(B,1)<gauss.nbgauss
fprintf('Nombre de points avant elimination  = %d\n',gauss.nbgauss);
gauss.coord=B;
w = gauss.w;
gauss.w = gauss.w(I);
K = setdiff(1:gauss.nbgauss,I);
for i=1:length(K)
gauss.w(J(K(i)))=gauss.w(J(K(i)))+w(K(i));    
end

gauss.nbgauss = length(gauss.w);
fprintf('Nombre de points apres elimination  = %d\n',gauss.nbgauss);
end
end
