function [Bsa,Ab,detJ] = calc_Bs(elem,xnode)
% function [Ba,Ab,detJ] = calc_Bs(elem,xnode)

nbelem = getnbelem(elem);
xnodeglobal = xnode;

gauss = calc_gauss(elem,0);
xgauss = gauss.coord;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
% DDNqlocal = getDDNq(elem,xgauss);

% SHEAR (UZ,RX,RY,A)=(UZ,-BY,BX,A)
Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(1,:)-xnode(3,:)];
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

C = Dxy(:,1)./l;
S = Dxy(:,2)./l;

t2 = [ Ji(1,1,:,:)^2 Ji(1,2,:,:)^2 2*Ji(1,1,:,:)*Ji(1,2,:,:)
       Ji(2,1,:,:)^2 Ji(2,2,:,:)^2 2*Ji(2,1,:,:)*Ji(2,2,:,:)
       Ji(1,1,:,:)*Ji(2,1,:,:) Ji(1,2,:,:)*Ji(2,2,:,:) Ji(1,1,:,:)*Ji(2,2,:,:)+Ji(1,2,:,:)*Ji(2,1,:,:) ];
T2 = [t2,zeros(3);zeros(3),t2];
% Ta = [DDNqlocal(:,4:6,:,:);DDNqlocal(:,4:6,:,:)].*[C,C,C,S,S,S]';
Ta = repmat([-8,0,0;0,0,-8;-4,4,-4],2,1).*[C,C,C,S,S,S]';
Pa = T2*Ta;

mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnodeglobal,gauss.coord);
Db = calc_opmatbending(mat,elem,xnodeglobal,gauss.coord);
Hb = [ Db(1,1) Db(3,3) 2*Db(1,3) Db(1,3) Db(2,3) Db(1,2)+Db(3,3)
       Db(1,3) Db(2,3) Db(1,2)+Db(3,3) Db(3,3) Db(2,2) 2*Db(2,3) ];

Bsa = Hb*Pa;

Aa = 2/3*l.*eye(3) - [Dxy(:,1),Dxy(:,2)]*(Hs\Bsa);        

a = -1;
b = 1/2*Dxy(:,1);
c = 1/2*Dxy(:,2);

Aw = zerosND([3,3*3,sizeND(detJ,1)]);
for i=1:3
    Aw(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a,b(i),c(i),-a,b(i),c(i)];
end

Ab = Aa\Aw;

Abtemp = Ab;
Ab(:,2:3:end) = -Abtemp(:,3:3:end); % RX = -BY
Ab(:,3:3:end) = Abtemp(:,2:3:end); % RY = BX
