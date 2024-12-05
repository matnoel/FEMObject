function bsd = righthandSD(GSD,b)

if ~israndom(b)
    bsd=b;
else
 tolfactSD = getparam(GSD,'finalSDfacttol');   
 
 if isa(b,'PCRADIALMATRIX')
     mm=length(b);
 else
     mm=min(length(getPC(b)),numel(b)); 
 end
 
 
bsd = spectral_decomposition(b,...
     'tol',getparam(GSD,'tol')*tolfactSD,...
     'nbfoncmax',mm);  
 fprintf('SD of right-hand side : %d -> %d\n',mm,length(bsd));
end 