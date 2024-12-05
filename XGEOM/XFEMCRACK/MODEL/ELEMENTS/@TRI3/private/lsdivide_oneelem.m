function [connecin,connecout,xnodeplus]=lsdivide_oneelem(ls,xnode)

if all(ls>=0)
connecin = zeros(0,3);
connecout = 1:3;
elseif all(ls<=0)
connecout = zeros(0,3);
connecin = 1:3;
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
    xnodetotal = [xnode;x4;x5];
    
    nisol = seg2node(1);
      
    % ------- 2 SEGMENTS CONTIGUS ------
    connectri = [nisol,5,4];
    polyconnec = [4,5,mod(nisol:nisol+1,3)+1];
    
    [subconnec,xcenter]=splitpoly(polyconnec,xnodetotal(polyconnec,:),6);
    
    if sign(ls(nisol))==1
     connecin = subconnec;
     connecout = connectri;
    else
     connecout = subconnec;
     connecin = connectri;
    end
    
    xnodeplus = [x4;x5;xcenter];
    
    
elseif length(segcut)==1 && length(nodecut)==1
% --------------------------------
% LA LEVELSET PASSE PAR 1 SEGMENT ET 1 SOMMET 
% --------------------------------
seg1node = segnode(:,segcut);
xseg1 = xnode(seg1node,:);
x4 = calcnodecut(xseg1(1,:),xseg1(2,:),ls(seg1node));
%nexclu = setdiff(1:4,[nodecut,seg1node']);
xnodetotal = [xnode;x4];

nisol = seg1node(1);

connectri1 = [nodecut,seg1node(1),4];
connectri2 = [nodecut,4,seg1node(2)];
xnodeplus = [x4];

    if sign(ls(seg1node(1)))==1
     connecin = connectri2;
     connecout = connectri1;
    else
     connecout =  connectri2;
     connecin = connectri1;
    end


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

