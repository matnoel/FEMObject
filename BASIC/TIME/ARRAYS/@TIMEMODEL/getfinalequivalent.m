function b0 = getfinalequivalent(T)
% function b0 = getfinalequivalent(T)

Mt = getMmatrix(T);
b0 = zeros(1,length(T));
b0(end) = 1 ; 
b0 = b0/Mt;
b0=TIMEMATRIX(b0,T,[1,1]);
