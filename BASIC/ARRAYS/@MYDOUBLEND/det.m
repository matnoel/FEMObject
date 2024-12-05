function J = det(J)
% function J = det(J)

detJ = det3D(reshape3D(J.double));
J.double = reshape(detJ,[1,1,sizeND(J.double)]);
