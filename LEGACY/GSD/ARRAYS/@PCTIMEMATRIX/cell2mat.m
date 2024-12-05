function u = cell2mat(u)

if isa(u.value,'cell')
    
nt = length(u);
s = size(u);
PC = getPC(u);
for k=1:nt
 u.value{k}=double(cell2mat(u.value{k}))';
end
 u.value = [u.value{:}]';
 u.value = PCMATRIX(u.value,[prod(s),nt],PC);

elseif isa(u.value,'PCRADIALMATRIX')
  u.value=cell2mat(u.value);
end