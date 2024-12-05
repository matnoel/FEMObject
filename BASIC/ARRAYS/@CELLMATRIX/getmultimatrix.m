function u=getmultimatrix(u,k)

if ~isa(k,'char') && ~strcmp(k,':') 
 u.value=u.value(:,k);
 u.sm = [length(k),1];
end