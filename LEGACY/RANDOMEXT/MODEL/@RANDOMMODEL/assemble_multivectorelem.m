function f = assemble_multivectorelem(S,fe,varargin)


liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if ~isa(fe,'cell') &  length(liste)>1
     error('specifier le groupe d''element') 
end

if isa(fe,'cell')
n=size(fe{liste(1)},2);
nmat=size(fe{liste(1)},4);
else
n=size(fe,2);
nmat=size(fe,4);    
end

f=sparse(S.nbddl,n*nmat);
for p=liste
  if ~isa(fe,'cell')
      fep = fe;
  else
      fep = fe{p};
  end
 [ip,valp,n,nmat]=multivectorelem(S.groupelem{p},double(fep)); 
 f(ip,:) = f(ip,:) + valp ; 
end

if nmat>1
f=reshape(f,S.nbddl*n,nmat);
f = MULTIMATRIX(f,[S.nbddl,n]);
end

