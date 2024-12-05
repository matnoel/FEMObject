function [B,detJ] = calc_B(elem,xnode,xgauss)
% function [B,detJ] = calc_B(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = getnbddl(elem);
xnodeglobal = xnode;

[detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xgauss);
N = MYDOUBLEND(getN(elem,xgauss));
DN = Ji*DNlocal;

% MEMBRANE (UX,UY)
Bm = zerosND([3,4*2,sizeND(DN)]);
Bm(1,1:2:end) = DN(1,:);
Bm(2,2:2:end) = DN(2,:);
Bm(3,2:2:end) = DN(1,:);
Bm(3,1:2:end) = DN(2,:);

% DRILLING RZ
Bd = zerosND([1,4,sizeND(N)]);
Bd(1,:) = N(1,:);

% BENDING + SHEAR (UZ,RX,RY)
Bb = zerosND([3,4*3,sizeND(DN)]);
Bb(1,3:3:end) = DN(1,:);
Bb(2,2:3:end) = -DN(2,:);
Bb(3,2:3:end) = -DN(1,:);
Bb(3,3:3:end) = DN(2,:);

Bs = zerosND([2,4*3,sizeND(DN)]);
Bs(1,1:3:end) = DN(1,:); % UZ'+RY
Bs(1,3:3:end) = N(1,:);
Bs(2,1:3:end) = DN(2,:); % UZ'-RX
Bs(2,2:3:end) = -N(1,:);

% TOTAL (UX,UY,UZ,RX,RY,RZ)
B = zerosND([9,nbddl,sizeND(detJ)]);
repm=[];repbs=[];repd=[];
for i=1:4
    repm = [repm,6*(i-1)+[1,2]];
    repbs = [repbs,6*(i-1)+[3,4,5]];
    repd = [repd,6*(i-1)+6];
end
B(1:3,repm) = Bm;
B(4:6,repbs) = Bb;
B(7:8,repbs) = Bs;
B(9,repd) = Bd;

% ROTATION
P = calc_P(elem);
B = B*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
