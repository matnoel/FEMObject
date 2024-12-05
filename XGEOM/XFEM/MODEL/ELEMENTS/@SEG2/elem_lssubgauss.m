function [subgaussin,subgaussout] = elem_lssubgauss(elem,ls,ordre,choix)
%function [subgaussin,subgaussout] = elem_lssubgauss(elem,ls,ordre,choix)
% si choix=2 : calcul des poids et coordonnees des sous-points de gauss
% si choix=1 : calcul des poids
if nargin==3
    choix = 2;
end


subgaussin.w=[];
subgaussin.coord=zeros(0,1);
subgaussout.w=[];
subgaussout.coord=zeros(0,1);

seg2gauss = gaussseg2(ordre);
if all(ls>=0)
    subgaussout=seg2gauss;
elseif all(ls<=0)
    subgaussin=seg2gauss; 
else
    
[conneclocalin,conneclocalout,xlnodeplus]=lsdivide_oneelem(ls);%

xlnodeplus = [nodelocalcoordseg2();xlnodeplus];

for k=1:size(conneclocalin,1)
xnode=xlnodeplus(conneclocalin(k,:),:);      
if choix==2
    c = Nseg2(seg2gauss.coord)*xnode;
    subgaussin.coord = [subgaussin.coord;c] ;
end
w = seg2gauss.w * detJseg2(xnode);
subgaussin.w = [subgaussin.w;w(:)] ;
end

subgaussin.nbgauss = numel(subgaussin.w);
subgaussout.nbgauss = numel(subgaussout.w);


for k=1:size(conneclocalout,1)
xnode=xlnodeplus(conneclocalout(k,:),:);      
if choix==2
c = Nseg2(seg2gauss.coord)*xnode;
subgaussout.coord = [subgaussout.coord;c] ;
end
w = seg2gauss.w * detJseg2(xnode);
subgaussout.w = [subgaussout.w;w(:)] ;
end


end
    