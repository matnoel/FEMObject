function varargout = ...
    elem_bilssubgauss(elem,ls1,ls2,ordre,choix)
%function [subgaussinin,subgaussinout,subgaussoutin,subgaussoutout] = elem_bilssubgauss(elem,ls1,ls2,ordre,choix)
% si choix=2 : calcul des poids et coordonnees des sous-points de gauss (par defaut)
% si choix=1 : calcul des poids
if nargin==4
    choix = 2;
end

subgauss = cell(1,4); 
for i=1:length(subgauss)
subgauss{i}.w = [];  
subgauss{i}.coord=zeros(0,2);
end

tri3gauss = gausstri3(ordre);

in1 = all(ls1<=0);
in2 = all(ls2<=0);
out1 = all(ls1>=0);
out2 = all(ls2>=0);

if out1 && out2
    subgauss{4}=tri3gauss;
elseif out1 && in2
    subgauss{3}=tri3gauss;
elseif in1 && out2
    subgauss{2}=tri3gauss;
elseif in1 && in2
    subgauss{1}=tri3gauss;
elseif in1
    [subgauss{1},subgauss{2}] = elem_lssubgauss(elem,ls2,ordre,choix);
elseif out1 
    [subgauss{3},subgauss{4}] = elem_lssubgauss(elem,ls2,ordre,choix);
elseif in2
    [subgauss{1},subgauss{3}] = elem_lssubgauss(elem,ls1,ordre,choix);
elseif out2
    [subgauss{2},subgauss{4}] = elem_lssubgauss(elem,ls1,ordre,choix);
else
subconnec=cell(1,4);
[subconnec{1},subconnec{2},subconnec{3},subconnec{4},xlnodeplus]=bilsdivide_oneelem(ls1,ls2);%
xlnodeplus = [nodelocalcoordtri3();xlnodeplus];

for i=1:length(subconnec)

for k=1:size(subconnec{i},1)
xnode=xlnodeplus(subconnec{i}(k,:),:); 
switch size(subconnec{i},2)
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
subgauss{i}.w = [subgauss{i}.w;w(:)] ;
if choix==2
subgauss{i}.coord = [subgauss{i}.coord;c] ;
end
end
subgauss{i}.nbgauss = numel(subgauss{i}.w);

end

end

varargout = subgauss;
    
function detJ = detJqua4(qua4ref,xnode,xi)

detJ=double(calc_detJ(qua4ref,xnode,xi));
detJ=detJ(:);

return
