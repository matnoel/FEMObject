function am=mean(apc)

Hm = mean(apc.POLYCHAOS);
if iscell(apc)    
rep = find(nonzeros(Hm));
am =  getmultimatrix(apc.MULTIMATRIX,rep); 
am =  getvalue(multimtimes(Hm(rep)',am));
am =  am{1};
else
am = reshape(double(apc)*Hm,size(apc));
end