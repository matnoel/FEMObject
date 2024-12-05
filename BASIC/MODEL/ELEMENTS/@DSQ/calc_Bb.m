function [B,detJ] = calc_Bb(elem,xnode,xgauss)
% function [B,detJ] = calc_Bb(elem,xnode,xgauss)

nbelem = getnbelem(elem);

[detJ,J,Ji,xnode,DNllocal] = calc_detJ(elem,xnode,xgauss);
DNqlocal = getDNq(elem,xgauss);
DNl = Ji*DNllocal;
DNq = Ji*DNqlocal;

% BENDING (UZ,RX,RY)=(UZ,-BY,BX)
Bb = zerosND([3,4*3,sizeND(detJ)]);
Bb(1,2:3:end) = DNl(1,:);
Bb(2,3:3:end) = DNl(2,:);
Bb(3,2:3:end) = DNl(2,:);
Bb(3,3:3:end) = DNl(1,:);

Bbtemp = Bb;
Bb(:,2:3:end) = -Bbtemp(:,3:3:end); % RX = -BY
Bb(:,3:3:end) = Bbtemp(:,2:3:end); % RY = BX

% BENDING A
Dxy = [xnode(2,:)-xnode(1,:);xnode(3,:)-xnode(2,:);xnode(4,:)-xnode(3,:);xnode(1,:)-xnode(4,:)];
l2 = (Dxy(:,1).^2 + Dxy(:,2).^2);
l = sqrt(l2);

C = Dxy(:,1)./l;
S = Dxy(:,2)./l;

DNZx = DNq(:,5:8).*C';
DNZy = DNq(:,5:8).*S';

Ba = zerosND([3,4,sizeND(detJ)]);
Ba(1,:) = DNZx(1,:);
Ba(2,:) = DNZy(2,:);
Ba(3,:) = DNZx(2,:)+DNZy(1,:);

% TOTAL (UZ,RX,RY,A)
B = zerosND([3,4*4,sizeND(detJ)]);
B = [Bb,Ba];
