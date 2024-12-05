function xnodeseg=contour_oneelem(ls,xnode)


if all(ls>=0) || all(ls<=0)
error('pas d''intersection de la levelset et de l''element')

elseif sum(ls==0)==2
    error('segment sur le bord, cas non traité');
    
else
    
segnode = [1,2,3;2,3,1];

segstate = sign(ls(segnode(1,:)).*ls(segnode(2,:)));
segcut = find(segstate==-1);
nodestate = sign(ls);
nodecut = find(nodestate==0);

if nargin==1
xnode=nodelocalcoordtri3();
end

if length(segcut)==2 
% --------------------------------
% LA LEVELSET PASSE PAR 2 SEGMENTS 
% --------------------------------
    if segcut(1)==1 && segcut(2)==3
        segcut=[3,1];
    end
    seg1node = segnode(:,segcut(1));
    seg2node = segnode(:,segcut(2));
    xseg1 = xnode(seg1node,:);
    x4 = calcnodecut(xseg1(1,:),xseg1(2,:),ls(seg1node));
    xseg2 = xnode(seg2node,:);
    x5 = calcnodecut(xseg2(1,:),xseg2(2,:),ls(seg2node));
    xnodeseg = [x4;x5];
    
    
elseif length(segcut)==1 && length(nodecut)==1
% --------------------------------
% LA LEVELSET PASSE PAR 1 SEGMENT ET 1 SOMMET 
% --------------------------------
seg1node = segnode(:,segcut);
xseg1 = xnode(seg1node,:);
x4 = calcnodecut(xseg1(1,:),xseg1(2,:),ls(seg1node));
%nexclu = setdiff(1:4,[nodecut,seg1node']);
xnodeseg = [xnode(nodecut,:);x4];

elseif length(nodecut)==2
xnodeseg = xnode(nodecut,:);    
    
elseif ~(isempty(segcut))
error('cas impossible')
else
error('pas prevu')

end


end

function x=calcnodecut(x1,x2,ls)
xi = ls(1)./(ls(1)-ls(2));
x = (1-xi)*x1 + xi*x2;
return

