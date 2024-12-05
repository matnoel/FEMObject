function [N,detJ,x] = calc_N(elem,xnode,xgauss)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss)

nbelem = getnbelem(elem);
nbddl = getnbddl(elem);
param = getparam(elem);
L = param.L;

Nu = getN(elem,xgauss);
Nv = getNv(elem,xgauss);
DNvlocal = getDNv(elem,xgauss);
[detJ,J,Ji] = calc_detJ(elem,xnode,xgauss);
DNv = Ji*DNvlocal;

switch getindim(elem)
    case 2
        repu = [1,4];
        repv = [2,3,5,6];
        N = zerosND([3,6,sizeND(Nv)]);
        N(1,repu) = Nu;
        N(2,repv) = Nv;
        N(3,repv) = DNv;
    case 3
        repux = [1,7];
        repvy = [2,6,8,12];
        repvz = [3,5,9,11];
        reptx = [4,10];
        N = zerosND([6,12,sizeND(Nv)]);
        T = diag([1,-1,1,-1]);
        Nvz = Nv*T ;
        DNvz = DNv*T ;
        N(1,repux) = Nu;
        N(2,repvy) = Nv;
        N(3,repvz) = Nvz;
        N(4,reptx) = Nu ;
        N(5,repvz) = -DNvz; %  RY = -UZ'
        N(6,repvy) = DNv; % RZ = UY'
end

% ROTATION
P = calc_P(elem);
N = N*P;

if nargout>2
    x = calc_x(elem,xnode,xgauss);
end
