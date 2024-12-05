function M = assemble_mixedmatrixelem(S1,S2,me,varargin)
% s taille de la matrice finale
liste = getcharin('selgroup',varargin,1:S1.nbgroupelem);
s=[S1.nbddl,S2.nbddl];
i=[];j=[];val = [];
for p=liste
  [ip,jp,valp]=mixedmatrixelem(S1.groupelem{p},S2.groupelem{p},double(me{p}),s); 
   i=[i;ip];j=[j;jp];val=[val;valp];
end

M=sparse(i,j,val,s(1),s(2));


