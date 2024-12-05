function [B,detJ,x] = calc_B(elem,xnode,xgauss)
% function [B,detJ,x] = calc_B(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = getnbddl(elem);

DNlocal = getDN(elem,xgauss);
DDNvlocal = getDDNv(elem,xgauss);

[detJ,J,Ji] = calc_detJ(elem,xnode,xgauss);

DN = Ji*DNlocal;
DDNv = Ji*Ji*DDNvlocal;

switch getindim(elem)
    case 2
        repu = [1,4];
        repv = [2,3,5,6];
        B = zerosND([2,nbddl,sizeND(DDNv)]);
        B(1,repu) = DN;
        B(2,repv) = DDNv;
    case 3
        repux = [1,7];
        repvy = [2,6,8,12];
        repvz = [3,5,9,11];
        reptx = [4,10];
        B = zerosND([4,nbddl,sizeND(DDNv)]);
        B(1,repux) = DN;
        B(2,reptx) = DN;
        B(4,repvy) = DDNv;
        DDNvz = DDNv;
        DDNvz(1,[2,4]) = -DDNvz(1,[2,4]); % RY = -UZ'
        B(3,repvz) = DDNvz;
end

% ROTATION
P = calc_P(elem);
B = B*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
