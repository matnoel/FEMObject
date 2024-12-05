function M = assemble_multimatrixelem(S,me,varargin)

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if ~isa(me,'cell') &  length(liste)>1
     error('specifier le groupe d''element') 
end

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


