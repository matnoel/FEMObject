function [subgaussin,subgaussout] = elem_lssubgauss(elem,ls,ordre,choix)
%function [subgaussin,subgaussout] = lssubgauss(elem,ls,ordre,choix)
% si choix=2 : calcul des poids et coordonnees des sous-points de gauss (par defaut)
% si choix=1 : calcul des poids
if nargin==3
    choix = 2;
end
subgaussin.w=[];
subgaussin.coord=zeros(0,2);
subgaussout.w=[];
subgaussout.coord=zeros(0,2);

tri3gauss = gausstri3(ordre);

if all(ls>=0)
    subgaussout=tri3gauss;
elseif all(ls<=0)
    subgaussin=tri3gauss;   
else
[conneclocalin,conneclocalout,xlnodeplus]=lsdivide_oneelem(ls);%

%tri3ref = TRI3();
%qua4ref = QUA4();
xlnodeplus = [nodelocalcoordtri3();xlnodeplus];

%elem_gauss(tri3ref,ordre);
%qua4gauss = elem_gauss(qua4ref,ordre);



for k=1:size(conneclocalin,1)
xnode=xlnodeplus(conneclocalin(k,:),:); 
switch size(conneclocalin,2)
    case 3       
if choix==2
    c = Ntri3(tri3gauss.coord)*xnode;
end
w = tri3gauss.w * detJtri3(xnode);
    case 4       
if choix==2
c = calc_x(qua4ref,xnode,qua4gauss.coord);
end
w = qua4gauss.w * detJqua4(xnode);
end
subgaussin.w = [subgaussin.w;w(:)] ;
if choix==2
subgaussin.coord = [subgaussin.coord;c] ;
end
end


for k=1:size(conneclocalout,1)
xnode=xlnodeplus(conneclocalout(k,:),:); 
switch size(conneclocalout,2)
    case 3       
if choix==2
c = Ntri3(tri3gauss.coord)*xnode;
end
w = tri3gauss.w * detJtri3(xnode);
    case 4       
if choix==2
c = calc_x(qua4ref,xnode,qua4gauss.coord);
end
w = qua4gauss.w * detJqua4(xnode);
end
subgaussout.w = [subgaussout.w;w(:)] ;
if choix==2
subgaussout.coord = [subgaussout.coord;c] ;
end
end

end

subgaussin.nbgauss = numel(subgaussin.w);
subgaussout.nbgauss = numel(subgaussout.w);

    
function detJ = detJqua4(qua4ref,xnode,xi)

detJ=double(calc_detJ(qua4ref,xnode,xi));
detJ=detJ(:);

return
