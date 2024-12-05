function [N,detJ,x] = calc_N(elem,xnode,xgauss)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = 3*6;
xnodeglobal = xnode;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
Nl = MYDOUBLEND(getN(elem,xgauss));
Nq = MYDOUBLEND(getNq(elem,xgauss));

Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(1,:)-xnode(3,:)];

% MEMBRANE (UX,UY,RZ)
Zmu = zerosND([3,3*3,sizeND(detJ)]);
Zmv = zerosND([3,3*3,sizeND(detJ)]);
for i=1:3
    Zmu(i,[3*(i-1)+1,mod(3*i,9)+1]) = 1/2;
    Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
    Zmu(i,mod(3*i,9)+3) = Dxy(i,2)./4;
    Zmv(i,[3*(i-1)+2,mod(3*i,9)+2]) = 1/2;
    Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
    Zmv(i,mod(3*i,9)+3) = -Dxy(i,1)./4;
end

NZmu = Nq(:,4:6)*Zmu;
NZmv = Nq(:,4:6)*Zmv;

Nm = zerosND([3,3*3,sizeND(detJ)]);
Nm(1,1:3:end) = Nq(1,1:3);
Nm(2,2:3:end) = Nq(1,1:3);
Nm(1,:) = Nm(1,:)+NZmu(1,:);
Nm(2,:) = Nm(2,:)+NZmv(1,:);
Nm(3,3:3:end) = Nl;

% Zmu = zerosND([3,3*3,sizeND(detJ)]);
% Zmv = zerosND([3,3*3,sizeND(detJ)]);
% for i=1:3
%     Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
%     Zmu(i,mod(3*i,9)+3) = Dxy(i,2)./4;
%     Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
%     Zmv(i,mod(3*i,9)+3) = -Dxy(i,1)./4;
% end
% 
% NZmu = Nq(:,4:6)*Zmu;
% NZmv = Nq(:,4:6)*Zmv;
% 
% Nm = zerosND([3,3*3,sizeND(detJ)]);
% Nm(1,1:3:end) = Nl;
% Nm(2,2:3:end) = Nl;
% Nm(1,:) = Nm(1,:)+NZmu(1,:);
% Nm(2,:) = Nm(2,:)+NZmv(1,:);
% Nm(3,3:3:end) = Nl;

% BENDING + SHEAR (UZ,RX,RY)=(UZ,-BY,BX)
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
Hs = calc_opmatshear(mat,elem,xnodeglobal,xgauss);
Db = calc_opmatbending(mat,elem,xnodeglobal,xgauss);
Hb = [ Db(1,1) Db(3,3) 2*Db(1,3) Db(1,3) Db(2,3) Db(1,2)+Db(3,3)
       Db(1,3) Db(2,3) Db(1,2)+Db(3,3) Db(3,3) Db(2,2) 2*Db(2,3) ];

Ba = Hb*Pa;

Aa = 2/3*l.*eye(3) - [Dxy(:,1),Dxy(:,2)]*(Hs\Ba);        

a = -1;
b = 1/2*Dxy(:,1);
c = 1/2*Dxy(:,2);

Aw = zerosND([3,3*3,sizeND(detJ)]);
for i=1:3
    Aw(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a,b(i),c(i),-a,b(i),c(i)];
end

Ab = Aa\Aw;

Zbsx = C.*Ab;
Zbsy = S.*Ab;

NZbsx = Nq(:,4:6)*Zbsx;
NZbsy = Nq(:,4:6)*Zbsy;

Nbs = zerosND([3,3*3,sizeND(detJ)]);
Nbs(1,1:3:end) = Nl;
Nbs(2,2:3:end) = Nl;
Nbs(3,3:3:end) = Nl;
Nbs(2,:) = Nbs(2,:)+NZbsx(1,:);
Nbs(3,:) = Nbs(3,:)+NZbsy(1,:);

Nbstemp = Nbs;
Nbs(:,2:3:end) = -Nbstemp(:,3:3:end); % RX = -BY
Nbs(:,3:3:end) = Nbstemp(:,2:3:end); % RY = BX

Nbstemp = Nbs;
Nbs(2,:) = -Nbstemp(3,:); % RX = -BY
Nbs(3,:) = Nbstemp(2,:); % RY = BX

% TOTAL (UX,UY,UZ,RX,RY,RZ)
N = zerosND([6,3*6,sizeND(detJ)]);
repm=[];repbs=[];
for i=1:3
    repm = [repm,6*(i-1)+[1,2,6]];
    repbs = [repbs,6*(i-1)+[3,4,5]];
end
N([1,2,6],repm) = Nm;
N([3,4,5],repbs) = Nbs;

% ROTATION
P = calc_P(elem);
N = N*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
