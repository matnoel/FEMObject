function [B,detJ]=calc_Blsout(elem,xnode,xgauss,ls)

lsxnode = double(ls);
lsxnode = lsxnode(getconnec(elem)');
lsxnode = reshape(lsxnode,[size(lsxnode,1),1,size(lsxnode,2)]);
Nuni = getN(elem,xgauss);
lsx = Nuni*(abs(lsxnode)-lsxnode);

nbelem=getnbelem(elem);
nbddl = getnbddl(elem);
n=getnbddlpernode(elem);
ng=getnbddlpergauss(elem);

[DN,detJ] = calc_DN(elem,xnode,xgauss);
dlsx = DN*(abs(lsxnode)-lsxnode);

switch getindim(elem)
    case 1
B = DN ;
    case 2
B=zerosND(3,nbddl,nbelem);
B(1,1:2:end,:)=DN(1,:,:);
B(2,2:2:end,:)=DN(2,:,:);
B(3,1:2:end,:)=DN(2,:,:);
B(3,2:2:end,:)=DN(1,:,:);
        dlsx = [dlsx(1),0;0,dlsx(2);dlsx(2),dlsx(1)];    
    case 3
B=zerosND(6,nbddl,nbelem);
B(1,1:3:end)=DN(1,:);
B(2,2:3:end)=DN(2,:);
B(3,3:3:end)=DN(3,:);
B(4,1:3:end)=DN(2,:);
B(4,2:3:end)=DN(1,:);
B(5,1:3:end)=DN(3,:);
B(5,3:3:end)=DN(1,:);
B(6,2:3:end)=DN(3,:);
B(6,3:3:end)=DN(2,:);
dlsx = [dlsx(1),0,0;0,dlsx(2),0;0,0,dlsx(3);dlsx(2),dlsx(1),0;dlsx(3),0,dlsx(1);0,dlsx(3),dlsx(2)];            
end


N = calc_DN(elem,xnode,xgauss);

B = [B ,lsx*B + dlsx*N];

Bout = B ; 
Bout(:,1:2:end) = B(:,1:nbddl/2);
Bout(:,2:2:end) = B(:,nbddl/2+1:end);
B=Bout;