function [B,detJ] = calc_B(elem,xnode,xgauss)
% function [B,detJ] = calc_B(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = getnbddl(elem);
xnodeglobal = xnode;

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
DNqlocal = getDNq(elem,xgauss);
Nl = getN(elem,xgauss);
DNl = Ji*DNllocal;
DNq = Ji*DNqlocal;

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

DNZmu = DNq(:,5:8)*Zmu;
DNZmv = DNq(:,5:8)*Zmv;

Bm = zerosND([3,4*3,sizeND(detJ)]);
Bm(1,1:3:end) = DNq(1,1:4);
Bm(2,2:3:end) = DNq(2,1:4);
Bm(3,1:3:end) = DNq(2,1:4);
Bm(3,2:3:end) = DNq(1,1:4);
Bm(1,:) = Bm(1,:)+DNZmu(1,:);
Bm(2,:) = Bm(2,:)+DNZmv(2,:);
Bm(3,:) = Bm(3,:)+DNZmu(2,:)+DNZmv(1,:);

% Zmu = zerosND([4,4*3,sizeND(detJ)]);
% Zmv = zerosND([4,4*3,sizeND(detJ)]);
% for i=1:4
%     Zmu(i,3*(i-1)+3) = -Dxy(i,2)./4;
%     Zmu(i,mod(3*i,12)+3) = Dxy(i,2)./4;
%     Zmv(i,3*(i-1)+3) = Dxy(i,1)./4;
%     Zmv(i,mod(3*i,12)+3) = -Dxy(i,1)./4;
% end
% 
% DNZmu = DNq(:,5:8)*Zmu;
% DNZmv = DNq(:,5:8)*Zmv;
% 
% Bm = zerosND([3,4*3,sizeND(detJ)]);
% Bm(1,1:3:end) = DNl(1,:);
% Bm(2,2:3:end) = DNl(2,:);
% Bm(3,1:3:end) = DNl(2,:);
% Bm(3,2:3:end) = DNl(1,:);
% Bm(1,:) = Bm(1,:)+DNZmu(1,:);
% Bm(2,:) = Bm(2,:)+DNZmv(2,:);
% Bm(3,:) = Bm(3,:)+DNZmu(2,:)+DNZmv(1,:);

% DRILLING (UX,UY,RZ)
Bd = zerosND([1,4*3,sizeND(detJ)]);
Bd(1,2:3:end) = DNq(1,1:4)./2;
Bd(1,1:3:end) = -DNq(2,1:4)./2;
Bd = Bd + DNZmv(1,:)./2;
Bd = Bd - DNZmu(2,:)./2;
Bd(1,3:3:end) = - Nl;

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

t2 = [ Ji_seg(1,1,:,:)^2 Ji_seg(1,2,:,:)^2 2*Ji_seg(1,1,:,:)*Ji_seg(1,2,:,:)
       Ji_seg(2,1,:,:)^2 Ji_seg(2,2,:,:)^2 2*Ji_seg(2,1,:,:)*Ji_seg(2,2,:,:)
       Ji_seg(1,1,:,:)*Ji_seg(2,1,:,:) Ji_seg(1,2,:,:)*Ji_seg(2,2,:,:) Ji_seg(1,1,:,:)*Ji_seg(2,2,:,:)+Ji_seg(1,2,:,:)*Ji_seg(2,1,:,:) ];
T2 = [t2,zeros(3);zeros(3),t2];
Ta = [DDNqlocal_seg(:,5:8,:,:);DDNqlocal_seg(:,5:8,:,:)].*[C,C,C,S,S,S]';
Pa = T2*Ta;

mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnodeglobal,xgauss_seg);
Db = calc_opmatbending(mat,elem,xnodeglobal,xgauss_seg);
Hb = [ Db(1,1) Db(3,3) 2*Db(1,3) Db(1,3) Db(2,3) Db(1,2)+Db(3,3)
       Db(1,3) Db(2,3) Db(1,2)+Db(3,3) Db(3,3) Db(2,2) 2*Db(2,3) ];

Bsb = Hb*Pb;
Bsa = Hb*Pa;

Bsb = sum(wgauss_seg*Bsb,4);
Bsa = sum(wgauss_seg*Bsa,4);

Aa = 2/3*l.*eye(4) - [Dxy(:,1),Dxy(:,2)]*(Hs\Bsa);

a = -1;
b = 1/2*Dxy(:,1);
c = 1/2*Dxy(:,2);

Aw = zerosND([4,4*3,sizeND(detJ)]);
for i=1:4
    Aw(i,[3*(i-1)+[1:3],mod(3*i,12)+[1:3]]) = -[a,b(i),c(i),-a,b(i),c(i)];
end
Aw = Aw + [Dxy(:,1),Dxy(:,2)]*(Hs\Bsb);

Ab = Aa\Aw;

Zbx = C.*Ab;
Zby = S.*Ab;

DNZbx = DNq(:,5:8)*Zbx;
DNZby = DNq(:,5:8)*Zby;

Bb = zerosND([3,4*3,sizeND(detJ)]);
Bb(1,2:3:end) = DNl(1,:);
Bb(2,3:3:end) = DNl(2,:);
Bb(3,2:3:end) = DNl(2,:);
Bb(3,3:3:end) = DNl(1,:);
Bb(1,:) = Bb(1,:)+DNZbx(1,:);
Bb(2,:) = Bb(2,:)+DNZby(2,:);
Bb(3,:) = Bb(3,:)+DNZbx(2,:)+DNZby(1,:);

Bbtemp = Bb;
Bb(:,2:3:end) = -Bbtemp(:,3:3:end); % RX = -BY
Bb(:,3:3:end) = Bbtemp(:,2:3:end); % RY = BX

Bs = zerosND([2,4*3,sizeND(detJ)]);
Bs = Hs\(Bsb + Bsa*Ab);

Bstemp = Bs;
Bs(:,2:3:end) = -Bstemp(:,3:3:end); % RX = -BY
Bs(:,3:3:end) = Bstemp(:,2:3:end); % RY = BX

% TOTAL (UX,UY,UZ,RX,RY,RZ)
B = zerosND([9,nbddl,sizeND(detJ)]);
repm=[];repbs=[];
for i=1:4
    repm = [repm,6*(i-1)+[1,2,6]];
    repbs = [repbs,6*(i-1)+[3,4,5]];
end
B(1:3,repm) = Bm;
B(4:6,repbs) = Bb;
B(7:8,repbs) = Bs;
B(9,repm) = Bd;

% ROTATION
P = calc_P(elem);
B = B*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
