function [N,detJ,x] = calc_N(elem,xnode,xgauss)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = 3*6;

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

% BENDING (UZ,RX,RY)=(UZ,-BY,BX)
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

a = -3/2*Dxy(:,1)./l2;
b = (1/4*Dxy(:,1).^2-1/2*Dxy(:,2).^2)./l2;
c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
d = -3/2*Dxy(:,2)./l2;
e = (1/4*Dxy(:,2).^2-1/2*Dxy(:,1).^2)./l2;

Zbx = zerosND([3,3*3,sizeND(detJ)]);
Zby = zerosND([3,3*3,sizeND(detJ)]);
for i=1:3
    Zbx(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a(i),b(i),c(i),-a(i),b(i),c(i)];
    Zby(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[d(i),c(i),e(i),-d(i),c(i),e(i)];
end

NZbx = Nq(:,4:6)*Zbx;
NZby = Nq(:,4:6)*Zby;

Nb = zerosND([3,3*3,sizeND(detJ)]);
Nb(1,1:3:end) = Nl;
Nb(2,2:3:end) = Nq(1,1:3);
Nb(3,3:3:end) = Nq(1,1:3);
Nb(2,:) = Nb(2,:)+NZbx(1,:);
Nb(3,:) = Nb(3,:)+NZby(1,:);

% a = -3/2*Dxy(:,1)./l2;
% b = 3/4*Dxy(:,1).^2./l2;
% c = 3/4*(Dxy(:,1).*Dxy(:,2))./l2;
% d = -3/2*Dxy(:,2)./l2;
% e = 3/4*Dxy(:,2).^2./l2;
% 
% Zbx = zerosND([3,3*3,sizeND(detJ)]);
% Zby = zerosND([3,3*3,sizeND(detJ)]);
% for i=1:3
%     Zbx(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[a(i),b(i),c(i),-a(i),b(i),c(i)];
%     Zby(i,[3*(i-1)+[1:3],mod(3*i,9)+[1:3]]) = -[d(i),c(i),e(i),-d(i),c(i),e(i)];
% end
% 
% NZbx = Nq(:,4:6)*Zbx;
% NZby = Nq(:,4:6)*Zby;
% 
% Nb = zerosND([3,3*3,sizeND(detJ)]);
% Nb(1,1:3:end) = Nl;
% Nb(2,2:3:end) = Nl;
% Nb(3,3:3:end) = Nl;
% Nb(2,:) = Nb(2,:)+NZbx(1,:);
% Nb(3,:) = Nb(3,:)+NZby(1,:);

Nbtemp = Nb;
Nb(:,2:3:end) = -Nbtemp(:,3:3:end); % RX = -BY
Nb(:,3:3:end) = Nbtemp(:,2:3:end); % RY = BX

Nbtemp = Nb;
Nb(2,:) = -Nbtemp(3,:); % RX = -BY
Nb(3,:) = Nbtemp(2,:); % RY = BX

% TOTAL (UX,UY,UZ,RX,RY,RZ)
N = zerosND([6,3*6,sizeND(detJ)]);
repm=[];repb=[];
for i=1:3
    repm = [repm,6*(i-1)+[1,2,6]];
    repb = [repb,6*(i-1)+[3,4,5]];
end
N([1,2,6],repm) = Nm;
N([3,4,5],repb) = Nb;

% ROTATION
P = calc_P(elem);
N = N*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
