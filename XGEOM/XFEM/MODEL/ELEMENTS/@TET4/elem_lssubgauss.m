function [subgaussin,subgaussout] = elem_lssubgauss(elem,ls,ordre,choix)
%function [subgaussin,subgaussout] = lssubgauss(elem,ls,ordre,choix)
% si choix=2 : calcul des poids et coordonnees des sous-points de gauss (par defaut)
% si choix=1 : calcul des poids
if nargin==3
    choix = 2;
end
subgaussin.w=[];
subgaussin.coord=zeros(0,3);
subgaussout.w=[];
subgaussout.coord=zeros(0,3);

tet4gauss = gausstet4(ordre);

if all(ls>=0)
    subgaussout=tet4gauss;
elseif all(ls<=0)
    subgaussin=tet4gauss;   
else
[conneclocalin,conneclocalout,xlnodeplus]=lsdivide_oneelem(ls);%

xlnodeplus = [nodelocalcoordtet4();xlnodeplus];


for k=1:size(conneclocalin,1)
xnode=xlnodeplus(conneclocalin(k,:),:); 
      
if choix==2
    c = Ntet4(tet4gauss.coord)*xnode;
end
w = tet4gauss.w * detJtet4(xnode);

subgaussin.w = [subgaussin.w;w(:)] ;
if choix==2
subgaussin.coord = [subgaussin.coord;c] ;
end
end


for k=1:size(conneclocalout,1)
xnode=xlnodeplus(conneclocalout(k,:),:); 
     
if choix==2
c = Ntet4(tet4gauss.coord)*xnode;
end
w = tet4gauss.w * detJtet4(xnode);

subgaussout.w = [subgaussout.w;w(:)] ;
if choix==2
subgaussout.coord = [subgaussout.coord;c] ;
end
end

end

subgaussin.nbgauss = numel(subgaussin.w);
subgaussout.nbgauss = numel(subgaussout.w);

    