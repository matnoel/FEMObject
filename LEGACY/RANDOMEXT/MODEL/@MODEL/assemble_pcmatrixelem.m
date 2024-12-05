function M = assemble_pcmatrixelem(S,me,varargin)

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if ~isa(me,'cell') &&  length(liste)>1
     error('specifier le groupe d''element') 
end

for p=liste
  if ~isa(me,'cell')
      mep = me;
  else
      mep = me{p};
  end
  if isa(mep,'PCMYDOUBLEND')   
  mep = permutestodim(mep,4);
  L = getL(mep);
  
  Mp{p} = assemble_multimatrixelem(S,getV(mep),'selgroup',p,'cell');
  if isa(L,'PCMATRIX')  
      if isa(Mp{p},'double')
  Mp{p} = PCRADIALMATRIX({Mp{p}},size(Mp{p}),L);
      else
  Mp{p} = PCRADIALMATRIX(Mp{p},size(Mp{p}),L);
      end
  
  elseif isa(L,'POLYCHAOS')
  Mp{p} = PCMATRIX(Mp{p},size(Mp{p}),L);  
  end
  elseif isa(mep,'MYDOUBLEND')   
  Mp{p} = assemble_multimatrixelem(S,mep,'selgroup',p);    
  else
      error('mauvais argument')
  end
end


M = Mp{liste(1)};

for p=liste(2:end)
M = M+Mp{p};
end
