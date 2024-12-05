function [B,detJ,DN]=calc_DNvect(elem,xnode,xgauss)

if getparam(elem,'initializeBN')
    B = getparam(elem,'B');
    detJ = getparam(elem,'detJ');
    return
end

nbelem=getnbelem(elem);

[DN,detJ] = calc_DN(elem,xnode,xgauss);

if getnbddlpernode(elem)==1
nbddl = getnbnode(elem);
B=DN;
elseif getnbddlpernode(elem)>=getindim(elem)
nbddl = getnbnode(elem)*getindim(elem);
switch getindim(elem)
    case 1
B=DN ; 
    case 2
B=zerosND([4,nbddl,sizeND(DN)]);
B(1,1:2:end)=DN(1,:);
B(2,2:2:end)=DN(2,:);
B(3,1:2:end)=DN(2,:);
B(4,2:2:end)=DN(1,:);
        
    case 3
error('pas fait')
B=zerosND([6,nbddl,sizeND(DN)]);
B(1,1:3:end)=DN(1,:);
B(2,2:3:end)=DN(2,:);
B(3,3:3:end)=DN(3,:);
B(4,1:3:end)=DN(2,:);
B(4,2:3:end)=DN(1,:);
B(5,1:3:end)=DN(3,:);
B(5,3:3:end)=DN(1,:);
B(6,2:3:end)=DN(3,:);
B(6,3:3:end)=DN(2,:);

        
    end
else
  error('calc_B non defini pour ce typ de ddl aux noeuds')  
end


