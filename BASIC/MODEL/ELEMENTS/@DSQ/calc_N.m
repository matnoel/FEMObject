function [N,detJ,x] = calc_N(elem,xnode,xgauss)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = 4*6;
xnodeglobal = xnode;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
Nl = MYDOUBLEND(getN(elem,xgauss));
Nq = MYDOUBLEND(getNq(elem,xgauss));

Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(4,:)-xnode(3,:);xnode(1,:)-xnode(4,:)];

% MEMBRANE (UX,UY,RZ)
Zmu = zerosND([4,4*3,sizeND(detJ)]);
Zmv = zerosND([4,4*3,sizeND(detJ)]);
for i=1:4
    Zmu(i,[3*(i-1)+1,mod(3*i,12)+1]) = 1/2;
    Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
    Zmu(i,mod(3*i,12)+3) = Dxy(i,2)./4;
    Zmv(i,[3*(i-1)+2,mod(3*i,12)+2]) = 1/2;
    Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
    Zmv(i,mod(3*i,12)+3) = -Dxy(i,1)./4;
end

NZmu = Nq(:,5:8)*Zmu;
NZmv = Nq(:,5:8)*Zmv;

Nm = zerosND([3,4*3,sizeND(detJ)]);
Nm(1,1:3:end) = Nq(1,1:4);
Nm(2,2:3:end) = Nq(1,1:4);
Nm(1,:) = Nm(1,:)+NZmu(1,:);
Nm(2,:) = Nm(2,:)+NZmv(1,:);
Nm(3,3:3:end) = Nl;

% Zmu = zerosND([4,4*3,sizeND(detJ)]);
% Zmv = zerosND([4,4*3,sizeND(detJ)]);
% for i=1:4
%     Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
%     Zmu(i,mod(3*i,12)+3) = Dxy(i,2)./4;
%     Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
%     Zmv(i,mod(3*i,12)+3) = -Dxy(i,1)./4;
% end
% 
% NZmu = Nq(:,5:8)*Zmu;
% NZmv = Nq(:,5:8)*Zmv;
% 
% Nm = zerosND([3,4*3,sizeND(detJ)]);
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

gauss_seg = calc_gausspoints(POLYLEGENDRE(),2);
gauss_seg.coord = gauss_seg.coord*[1 1];
gauss_seg = permutegaussND(gauss_seg);
xgauss_seg = gauss_seg.coord;
wgauss_seg = gauss_seg.w;

[detJ_seg,J_seg,Ji_seg] = calc_detJ(elem,xnodeglobal,xgauss_seg);
DDNqlocal_seg = getDDNq(elem,xgauss_seg);

a = 1/2*Ji_seg(1,1,:,:)*Ji_seg(1,2,:,:);
b = 1/2*Ji_seg(2,1,:,:)*Ji_seg(2,2,:,:);
c = 1/4*(Ji_seg(1,1,:,:)*Ji_seg(2,2,:,:)+Ji_seg(1,2,:,:)*Ji_seg(2,1,:,:));

Pb = zerosND([2*3,4*3,sizeND(Ji_seg)]);
Pb(1:3,2:3:end) = [a;b;c]*[1,-1,1,-1];
Pb(4:6,3:3:end) = [a;b;c]*[1,-1,1,-1];
Pb = sum(wgauss_seg*Pb,4);

t2 = [ Ji_seg(1,1,:,:)^2 Ji_seg(1,2,:,:)^2 2*Ji_seg(1,1,:,:)*Ji_seg(1,2,:,:)
       Ji_seg(2,1,:,:)^2 Ji_seg(2,2,:,:)^2 2*Ji_seg(2,1,:,:)*Ji_seg(2,2,:,:)
       Ji_seg(1,1,:,:)*Ji_seg(2,1,:,:) Ji_seg(1,2,:,:)*Ji_seg(2,2,:,:) Ji_seg(1,1,:,:)*Ji_seg(2,2,:,:)+Ji_seg(1,2,:,:)*Ji_seg(2,1,:,:) ];
T2 = [t2,zeros(3);zeros(3),t2];
Ta = [DDNqlocal_seg(:,5:8,:,:);DDNqlocal_seg(:,5:8,:,:)].*[C,C,C,S,S,S]';
Pa = T2*Ta;
Pa = sum(wgauss_seg*Pa,4);

mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnodeglobal,gauss_seg.coord);
Db = calc_opmatbending(mat,elem,xnodeglobal,gauss_seg.coord);
Hb = [ Db(1,1) Db(3,3) 2*Db(1,3) Db(1,3) Db(2,3) Db(1,2)+Db(3,3)
       Db(1,3) Db(2,3) Db(1,2)+Db(3,3) Db(3,3) Db(2,2) 2*Db(2,3) ];

Bb = Hb*Pb;
Ba = Hb*Pa;

Aa = 2/3*l.*eye(4) - [Dxy(:,1),Dxy(:,2)]*(Hs\Ba);        

a = -1;
b = 1/2*Dxy(:,1);
c = 1/2*Dxy(:,2);

Aw = zerosND([4,4*3,sizeND(detJ)]);
for i=1:4
    Aw(i,[3*(i-1)+[1:3],mod(3*i,12)+[1:3]]) = -[a,b(i),c(i),-a,b(i),c(i)];
end
Aw = Aw + [Dxy(:,1),Dxy(:,2)]*(Hs\Bb);

Ab = Aa\Aw;

Zbsx = C.*Ab;
Zbsy = S.*Ab;

NZbsx = Nq(:,5:8)*Zbsx;
NZbsy = Nq(:,5:8)*Zbsy;

Nbs = zerosND([3,4*3,sizeND(detJ)]);
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
N = zerosND([6,4*6,sizeND(detJ)]);
repm=[];repbs=[];
for i=1:4
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
