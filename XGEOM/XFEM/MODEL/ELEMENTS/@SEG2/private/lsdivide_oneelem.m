function [connecin,connecout,xnodeplus]=lsdivide_oneelem(ls,xnode)

segcut = (sign(ls(1).*ls(2))==-1);
nodestate = sign(ls);
nodecut = find(nodestate==0);

if nargin==1
xnode=nodelocalcoordseg2();
end

if segcut
% --------------------------------
% LA LEVELSET PASSE PAR 2 SEGMENTS 
% --------------------------------
    x3 = calcnodecut(xnode(1,:),xnode(2,:),ls);
    xnodetotal = [xnode;x3];
    
    if sign(ls(1))==1
     connecin = [3,2];
     connecout = [1,3];
    else
     connecout = [3,2];
     connecin = [1,3];
    end
    
    xnodeplus = x3;
 
else
    connecin = [];
    connecout = [];
    if nargin==2
    xnodeplus = zeros(0,size(xnode,2));
    else
    xnodeplus = zeros(0,1);
    end

end


function x=calcnodecut(x1,x2,ls)
xi = ls(1)./(ls(1)-ls(2));
x = (1-xi)*x1 + xi*x2;
return
