function b0 = getinitialequivalent(T)
% function b0 = getinitialequivalent(T)

Mt = getMmatrix(T);
b0 = zeros(1,length(T));
b0(1) = 1 ; 
b0 = b0/Mt;
b0=TIMEMATRIX(b0,T,[1,1]);
