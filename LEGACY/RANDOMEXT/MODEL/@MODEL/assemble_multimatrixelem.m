function M = assemble_multimatrixelem(S,me,varargin)

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if ~isa(me,'cell') &  length(liste)>1
     error('specifier le groupe d''element') 
end

if ischarin('cell',varargin)
    
 for p=liste
  if ~isa(me,'cell')
      mep = me;
  else
      mep = me{p};
  end
  
  
  if isa(mep,'MYDOUBLEND') && size(mep,4)>1
  mep = MULTIMYDOUBLEND(mep,4);
  end
  
  if isa(mep,'MULTIMYDOUBLEND')  
  for k=1:numelm(mep)
  [ip,jp,valp]=matrixelem(S.groupelem{p},double(mep{k}),[S.nbddl,S.nbddl]); 
  Mk{k} = sparse(ip(:),jp(:),valp(:),S.nbddl,S.nbddl);   
  end
  Mp{p} = MULTIMATRIX(Mk,[S.nbddl,S.nbddl],[numelm(mep),1]);
  elseif isa(mep,'MYDOUBLEND')    
  [ip,jp,valp]=matrixelem(S.groupelem{p},double(mep),[S.nbddl,S.nbddl]); 
  Mp{p} =  sparse(ip(:),jp(:),valp(:),S.nbddl,S.nbddl);   
  end

 end
 
 
M = Mp{liste(1)};
for p=liste(2:end)
M = M+Mp{p};
end

    
else

i=[];j=[];val = [];
for p=liste
  if ~isa(me,'cell')
      mep = me;
  else
      mep = me{p};
  end
  if isa(mep,'MULTIMYDOUBLEND')
      mep = permutemultidim(mep,4);
  end
      
  [ip,jp,valp,nmat]=multimatrixelem(S.groupelem{p},double(mep),S); 
   i=[i;ip];j=[j;jp];val=[val;valp];
end

if nmat==1
M=sparse(i(:),j(:),val(:),S.nbddl,S.nbddl);
else
M=sparse(i(:),j(:),val(:),S.nbddl*S.nbddl,nmat);    
M = MULTIMATRIX(M,[S.nbddl,S.nbddl]);
end


end