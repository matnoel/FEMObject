function [connecin,connecout,xnodeplus]=lsdivide_oneelem(ls,xnode)

if all(ls>=0)
connecin = zeros(0,4);
connecout = 1:4;
elseif all(ls<=0)
connecout = zeros(0,4);
connecin = 1:4;
else
    
segnode = [1,2;2,3;1,3;1,4;2,4;3,4];

segstate = sign(ls(segnode(:,1)).*ls(segnode(:,2)));
segcut = find(segstate==-1);
nodestate = sign(ls);
nodecut = find(nodestate==0);

if nargin==1
xnode=nodelocalcoordtet4();
end

if length(segcut)==3
% --------------------------------
% LA LEVELSET PASSE PAR 3 SEGMENTS 
% --------------------------------
xnodeplus=zeros(3,3);
for k=1:3
segknode=segnode(segcut(k),:);
xsegk = xnode(segknode,:);
xnodeplus(k,:) = calcnodecut(xsegk(1,:),xsegk(2,:),ls(segknode));
end
xnodetotal = [xnode;xnodeplus];

[tet,polyconnec,n1]=polyhedracase3seg(segcut);
[subtet] = splitprisme(polyconnec,xnodetotal(polyconnec,:));
 
    if sign(ls(n1))==1
     connecin = subtet;
     connecout = tet;
    else
     connecout = subtet;
     connecin = tet;
    end
   
elseif length(segcut)==4
% --------------------------------
% LA LEVELSET PASSE PAR 4 SEGMENTS 
% --------------------------------    
xnodeplus=zeros(4,3);
for k=1:4
segknode=segnode(segcut(k),:);
xsegk = xnode(segknode,:);
xnodeplus(k,:) = calcnodecut(xsegk(1,:),xsegk(2,:),ls(segknode));
end
xnodetotal = [xnode;xnodeplus];
[polyconnec1,polyconnec2,n1]=polyhedracase4seg(segcut);
subtet1 = splitprisme(polyconnec1,xnodetotal(polyconnec1,:));    
subtet2 = splitprisme(polyconnec2,xnodetotal(polyconnec2,:));  
    if sign(ls(n1))==1
     connecin = subtet2;
     connecout = subtet1;
    else
     connecout = subtet2;
     connecin = subtet1;
    end
   
elseif length(segcut)==2 && length(nodecut)==1
% --------------------------------
% LA LEVELSET PASSE PAR 2 SEGMENTS 
% --------------------------------    
xnodeplus=zeros(2,3);
for k=1:2
segknode=segnode(segcut(k),:);
xsegk = xnode(segknode,:);
xnodeplus(k,:) = calcnodecut(xsegk(1,:),xsegk(2,:),ls(segknode));
end    
xnodetotal = [xnode;xnodeplus];

[tet,polyconnec,n1]=polyhedracase2seg(segcut,nodecut);
[subtet] = splitpyr(polyconnec);

    if sign(ls(n1))==1
     connecin = subtet;
     connecout = tet;
    else
     connecout = subtet;
     connecin = tet;
    end
    
elseif length(segcut)==1 && length(nodecut)==2
% --------------------------------
% LA LEVELSET PASSE PAR 1 SEGMENT 
% --------------------------------
xnodeplus=zeros(1,3);
for k=1:1
segknode=segnode(segcut(k),:);
xsegk = xnode(segknode,:);
xnodeplus(k,:) = calcnodecut(xsegk(1,:),xsegk(2,:),ls(segknode));
end    
xnodetotal = [xnode;xnodeplus];

[subtet1,subtet2,n1]=polyhedracase1seg(segcut);    
    if sign(ls(n1))==1
     connecin = subtet2;
     connecout = subtet1;
    else
     connecout = subtet2;
     connecin = subtet1;
    end    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% ADDITIONAL FUNS  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=calcnodecut(x1,x2,ls)
xi = ls(1)./(ls(1)-ls(2));
x = (1-xi)*x1 + xi*x2;
return


function [poly1,poly2,nisol]=polyhedracase3seg(segcut)

[temp,k] = ismember(segcut(:)',[1,2,5;2,3,6;4,5,6;1,3,4],'rows');
switch k
    case 1
       nisol=2;
       poly1 = [2,6,5,7];
       poly2 = [1,3,4,5,6,7];
    case 2
       nisol=3; 
       poly1 = [3,6,5,7];
       poly2 = [1,2,4,6,5,7];
    case 3
       nisol=4; 
       poly1 = [4,5,6,7];
       poly2 = [1,2,3,5,6,7];
    case 4
       nisol=1; 
       poly1 = [1,5,6,7];
       poly2 = [2,3,4,5,6,7];
    otherwise
        error('pas prevu')
end

return

function [poly1,poly2,n1]=polyhedracase4seg(segcut)

[temp,k] = ismember(segcut(:)',[1,2,4,6;1,3,5,6;2,3,4,5],'rows');
switch k
    case 1
       n1=2;
       poly1 = [2,4,6,8,7,5];
       poly2 = [3,1,8,7,5,6];
    case 2
       n1=1; 
       poly1 = [1,4,5,7,8,6];
       poly2 = [2,3,5,6,8,7];
    case 3
       n1=3; 
       poly1 = [3,4,6,7,8,5];
       poly2 = [1,2,6,5,8,7];
    otherwise
        error('pas prevu')
end

poly1 = poly1([1,6,3,2,5,4]);
poly2 = poly2([1,6,3,2,5,4]);

return

function [poly1,poly2,nisol]=polyhedracase2seg(segcut,nodecut)

switch nodecut
    case 1
     [temp,k] = ismember(segcut(:)',[2,6;2,5;5,6],'rows');
    
     switch k
        case 1
       nisol = 3;
       poly1=[1,5,6,3];
       poly2=[1,2,5,6,4];
       
        case 2
       nisol = 2;
       poly1=[1,5,6,2];
       poly2=[1,3,4,6,5];
        
        case 3
       nisol = 4;
       poly1=[1,5,6,4];
       poly2=[1,2,3,6,5];
            
     end
    case 2
      [temp,k] = ismember(segcut(:)',[3,4;3,6;4,6],'rows');
    
     switch k
        case 1
       nisol = 1;
       poly1=[2,5,6,1];
       poly2=[2,3,4,6,5];
       
        case 2
       nisol = 3;
       poly1=[2,5,6,3];
       poly2=[2,5,6,4,1];
        
        case 3
       nisol = 4;
       poly1=[2,6,4,5];
       poly2=[2,3,6,5,1];
            
     end
    case 3
       [temp,k] = ismember(segcut(:)',[1,4;1,5;4,5],'rows');
    
     switch k
        case 1
       nisol = 1;
       poly1=[3,1,5,6];
       poly2=[3,6,5,2,4];
       
        case 2
       nisol = 2;
       poly1=[3,6,5,2];
       poly2=[3,4,1,5,6];
        
        case 3
       nisol = 4;
       poly1=[3,4,5,6];
       poly2=[3,1,5,6,2];
            
     end
    case 4
       [temp,k] = ismember(segcut(:)',[1,2;1,3;2,3],'rows');
    
     switch k
        case 1
       nisol = 2;
       poly1=[4,5,2,6];
       poly2=[4,5,6,3,1];
       
        case 2
       nisol = 1;
       poly1=[4,1,5,6];
       poly2=[4,5,2,3,6];
        
        case 3
       nisol = 3;
       poly1=[4,5,3,6];
       poly2=[4,1,2,5,6];
            
     end
    otherwise
        error('pas prevu')
end

return


function [poly1,poly2,n1]=polyhedracase1seg(segcut)

switch segcut
    case 1
       n1=2;
       poly1 = [3,4,5,2];
       poly2 = [3,4,5,1];
    case 2
       n1=2; 
       poly1 = [1,4,5,2];
       poly2 = [1,4,5,3];
    case 3
       n1=1; 
       poly1 = [2,4,5,1];
       poly2 = [2,4,5,3];
    case 4
       n1=4; 
       poly1 = [2,3,5,4];
       poly2 = [2,3,5,1]; 
    case 5
       n1=4; 
       poly1 = [1,3,5,4];
       poly2 = [1,3,5,2];
    case 6
       n1=4; 
       poly1 = [1,2,5,4];
       poly2 = [1,2,5,3]; 
    otherwise
        error('pas prevu')
end


return

function tet = splitprisme(prisme,xnode)

   tet=[1,2,3,4;6,4,3,5;2,3,4,5];
   tet = prisme(tet);
return
   

function tet = splitpyr(pyr,xnode)

   tet=[1,2,3,5;1,3,4,5];
   tet = pyr(tet);
 
return
       
    