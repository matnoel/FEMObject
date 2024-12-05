function [connecin,connecout,xnodeplus]=lsdivide_oneelem(ls,xnode)

if all(ls>=0)
connecin = zeros(0,4);
connecout = 1:4;
elseif all(ls<=0)
connecout = zeros(0,4);
connecin = 1:4;
else
    
segnode = [1,2,3,4;2,3,4,1];

segstate = sign(ls(segnode(1,:)).*ls(segnode(2,:)));
segcut = find(segstate==-1);
nodestate = sign(ls);
nodecut = find(nodestate==0);

if nargin==1
xnode=nodelocalcoordqua4();;
end

if length(segcut)==2 
% --------------------------------
% LA LEVELSET PASSE PAR 2 SEGMENTS 
% --------------------------------
    if segcut(1)==1 & segcut(2)==4
        segcut=[4,1];
    end
    seg1node = segnode(:,segcut(1));
    seg2node = segnode(:,segcut(2));
    xseg1 = xnode(seg1node,:);
    x5 = calcnodecut(xseg1(1,:),xseg1(2,:),ls(seg1node));
    xseg2 = xnode(seg2node,:);
    x6 = calcnodecut(xseg2(1,:),xseg2(2,:),ls(seg2node));
    xnodetotal = [xnode;x5;x6];
    if mod(segcut(2)-segcut(1),4)==1
    nisol = seg2node(1);
    else
    nisol=[];    
    end
%    nisol = intersect(seg1node,seg2node);
    
    if ~isempty(nisol)
    % ------- 2 SEGMENTS CONTIGUS ------
    connectri = [nisol,6,5];
    polyconnec = [5,6,mod(nisol:nisol+2,4)+1];
    
    [subconnec,xcenter]=splitpoly(polyconnec,xnodetotal(polyconnec,:),7);
    
    if sign(ls(nisol))==1
     connecin = subconnec;
     connecout = connectri;
    else
     connecout = subconnec;
     connecin = connectri;
    end
    
    xnodeplus = [x5;x6;xcenter];
    else
    % ------- 2 SEGMENTS OPPOSES  ------- 
    poly1 = [seg1node(1),5,6,seg2node(2)];
    poly2 = [seg1node(2),seg2node(1),6,5];
    [subconnec1,xc1]=splitpoly(poly1,xnodetotal(poly1,:),7);        
    [subconnec2,xc2]=splitpoly(poly2,xnodetotal(poly2,:),8);        
    if sign(ls(seg1node(1)))==1
     connecin = subconnec2;
     connecout = subconnec1;
    else
     connecin = subconnec1;
     connecout = subconnec2;
    end
    
    xnodeplus = [x5;x6;xc1;xc2];
    end
    
elseif length(segcut)==1 & length(nodecut)==1
% --------------------------------
% LA LEVELSET PASSE PAR 1 SEGMENT ET 1 SOMMET 
% --------------------------------
seg1node = segnode(:,segcut);
xseg1 = xnode(seg1node,:);
x5 = calcnodecut(xseg1(1,:),xseg1(2,:),ls(seg1node));
%nexclu = setdiff(1:4,[nodecut,seg1node']);
xnodetotal = [xnode;x5];
if (mod(nodecut,4)+1)==seg1node(1)
nisol = seg1node(1);
polyconnec = [nodecut,5,seg1node(2),mod(seg1node(2),4)+1];
connectri = [nodecut,nisol,5];
[subconnec,xcenter]=splitpoly(polyconnec,xnodetotal(polyconnec,:),6);
else
nisol = seg1node(2);
polyconnec = [seg1node(1),5,nodecut,mod(nodecut,4)+1];
connectri = [nodecut,5,nisol];
[subconnec,xcenter]=splitpoly(polyconnec,xnodetotal(polyconnec,:),6);
end
xnodeplus = [x5;xcenter];

    if sign(ls(nisol))==1
     connecin = subconnec;
     connecout = connectri;
    else
     connecout = subconnec;
     connecin = connectri;
    end


elseif length(nodecut)==2 & (nodecut(2)-nodecut(1))==2
% --------------------------------
% LA LEVELSET PASSE PAR 2 SOMMETS OPPOSES
% --------------------------------
xnodeplus = zeros(0,2);
subconnec1 = [mod(nodecut(1),4)+1,nodecut(2),nodecut(1)];
subconnec2 = [mod(nodecut(2),4)+1,nodecut(1),nodecut(2)];
   if sign(ls(mod(nodecut(1),4)+1))==1
     connecin = subconnec2;
     connecout = subconnec1;
    else
     connecin = subconnec1;
     connecout = subconnec2;
    end

elseif ~(length(segcut)==0 )
    error('cas impossible')
else
    error('pas prevu')

end

end

function x=calcnodecut(x1,x2,ls)
xi = ls(1)./(ls(1)-ls(2));
x = (1-xi)*x1 + xi*x2;
return

function [subconnec,xcenter]=splitpoly(polyconnec,xpoly,numcenter)
n = length(polyconnec);

xcenter = sum(xpoly,1)/n;
num1 = numcenter*ones(n,1);
num2 = polyconnec(1:n);
num3 = polyconnec([2:n,1]);
subconnec = [num1(:),num2(:),num3(:)];

return