function f = assemble_pcvectorelem(S,fe,varargin)

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if ~isa(fe,'cell') &  length(liste)>1
     error('specifier le groupe d''element') 
end

for p=liste
  if ~isa(fe,'cell')
      fep = fe;
  else
      fep = fe{p};
  end
  if isa(fep,'PCMYDOUBLEND')   
  fep = permutestodim(fep,4);
  L = getL(fep);
  fp{p} = assemble_multivectorelem(S,getV(fep),'selgroup',p);
  if isa(L,'PCMATRIX')      
  fp{p} = PCRADIALMATRIX(fp{p},size(fp{p}),L);
  elseif isa(L,'POLYCHAOS')
  fp{p} = PCMATRIX(fp{p},size(fp{p}),L);  
  end
  elseif isa(fep,'MYDOUBLEND')   
  fp{p} = assemble_multivectorelem(S,fep,'selgroup',p);    
  else
      error('mauvais argument')
  end
end


f = fp{liste(1)};
for p=liste(2:end)
f = f+fp{p};
end
