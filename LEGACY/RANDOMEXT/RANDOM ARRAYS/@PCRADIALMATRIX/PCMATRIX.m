function A = PCMATRIX(R)

PC = R.POLYCHAOS ;
s  = size(R.V);
L =  R.D * double(R.L);
if iscell(R.V)
%warning('transformation en MULTIMATRIX double')
R.V = cell2mat(R.V);
%A = multimtimes(L',R.V);
A = double(R.V) * L;    

else
A = double(R.V) * L;    
end
A = PCMATRIX(A,s,PC);