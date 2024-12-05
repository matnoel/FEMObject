function apc = pcexpand(a)


PCTP = a.POLYCHAOSTP;
PC = POLYCHAOS(PCTP);
PC = calc_indices(PC);

ind = getindices(PC);
phi = a.phi ; 

for k=1:getnbdim(a)
    if ~isranddim(a,k)
   phi{k} = phi{k}*mean(PCTP,k);  
    end
   phi{k} =  phi{k}(ind(:,k)+1);
   
end
phi = prod([phi{:}],2);
apc = a.phi0.*PCMATRIX(phi,[1,1],PC); 

