function Ji = inv(J)
% function Ji = inv(J)

Ji = J;

detJ = det3D(reshape3D(J));
tcomJ = tcom3D(reshape3D(J));
[tcomJ,detJ] = samesize2D(tcomJ,detJ);

Ji.double = tcomJ./detJ;
Ji.double = reshape(Ji,size(J));
