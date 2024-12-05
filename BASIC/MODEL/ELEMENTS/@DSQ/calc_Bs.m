function [Bsb,Bsa,Ab,detJ] = calc_Bs(elem,xnode)
% function [Bsb,Bsa,Ab,detJ] = calc_Bs(elem,xnode)

nbelem = getnbelem(elem);
xnodeglobal = xnode;

gauss = calc_gausspoints(POLYLEGENDRE(),2);
gauss.coord = gauss.coord*[1 1];
gauss = permutegaussND(gauss);
xgauss = gauss.coord;
wgauss = gauss.w;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
DDNqlocal = getDDNq(elem,xgauss);

% SHEAR (UZ,RX,RY,A)=(UZ,-BY,BX,A)
Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(4,:)-xnode(3,:);xnode(1,:)-xnode(4,:)];
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

C = Dxy(:,1)./l;
S = Dxy(:,2)./l;

a = 1/2*Ji(1,1,:,:)*Ji(1,2,:,:);
b = 1/2*Ji(2,1,:,:)*Ji(2,2,:,:);
c = 1/4*(Ji(1,1,:,:)*Ji(2,2,:,:)+Ji(1,2,:,:)*Ji(2,1,:,:));

Pb = zerosND([2*3,4*3,sizeND(Ji)]);
Pb(1:3,2:3:end) = [a;b;c]*[1,-1,1,-1];
Pb(4:6,3:3:end) = [a;b;c]*[1,-1,1,-1];

t2 = [ Ji(1,1,:,:)^2 Ji(1,2,:,:)^2 2*Ji(1,1,:,:)*Ji(1,2,:,:)
       Ji(2,1,:,:)^2 Ji(2,2,:,:)^2 2*Ji(2,1,:,:)*Ji(2,2,:,:)
       Ji(1,1,:,:)*Ji(2,1,:,:) Ji(1,2,:,:)*Ji(2,2,:,:) Ji(1,1,:,:)*Ji(2,2,:,:)+Ji(1,2,:,:)*Ji(2,1,:,:) ];
T2 = [t2,zeros(3);zeros(3),t2];
Ta = [DDNqlocal(:,5:8,:,:);DDNqlocal(:,5:8,:,:)].*[C,C,C,S,S,S]';
Pa = T2*Ta;

mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnodeglobal,xgauss);
Db = calc_opmatbending(mat,elem,xnodeglobal,xgauss);
Hb = [ Db(1,1) Db(3,3) 2*Db(1,3) Db(1,3) Db(2,3) Db(1,2)+Db(3,3)
       Db(1,3) Db(2,3) Db(1,2)+Db(3,3) Db(3,3) Db(2,2) 2*Db(2,3) ];

Bsb = Hb*Pb;
Bsa = Hb*Pa;

Bsb = sum(wgauss*Bsb,4);
Bsa = sum(wgauss*Bsa,4);

Aa = 2/3*l.*eye(4) - [Dxy(:,1),Dxy(:,2)]*(Hs\Bsa);        

a = -1;
b = 1/2*Dxy(:,1);
c = 1/2*Dxy(:,2);

Aw = zerosND([4,4*3,sizeND(detJ,1)]);
for i=1:4
    Aw(i,[3*(i-1)+[1:3],mod(3*i,12)+[1:3]]) = -[a,b(i),c(i),-a,b(i),c(i)];
end
Aw = Aw + [Dxy(:,1),Dxy(:,2)]*(Hs\Bsb);

Ab = Aa\Aw;

Abtemp = Ab;
Ab(:,2:3:end) = -Abtemp(:,3:3:end); % RX = -BY
Ab(:,3:3:end) = Abtemp(:,2:3:end); % RY = BX

Bbtemp = Bsb;
Bsb(:,2:3:end) = -Bbtemp(:,3:3:end); % RX = -BY
Bsb(:,3:3:end) = Bbtemp(:,2:3:end); % RY = BX
