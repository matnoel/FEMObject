function Z=expand(apcr)

if ~iscell(apcr.V) && ~iscell(apcr.L)
Z = double(apcr.V)* apcr.D * double(apcr.L) ;
Z = PCMATRIX(Z,size(apcr),apcr.POLYCHAOS);
else 
PC = getPC(apcr.L);
L = getmultimatrix(apcr.L);
L = apcr.D*L;
N = length(PC);
Z = cell(N,1);

for k=1:N 
Lk = L{k};
Z{k} = double(multimtimes(Lk',apcr.V));
end

Z = PCMATRIX(Z,size(apcr),PC);

end
