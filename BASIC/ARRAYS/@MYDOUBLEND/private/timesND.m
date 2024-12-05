function C = timesND(A,B)
% function C = timesND(A,B)

C = times3D(reshape3D(A),reshape3D(B));
C = reshape(C,[size(A,2),size(B,1),sizeND(A)]);
