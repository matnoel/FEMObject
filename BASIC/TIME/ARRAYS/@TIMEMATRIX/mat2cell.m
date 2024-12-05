function u=mat2cell(u)

if ~iscell(u.value)
nt = getnt(u);
v = cell(1,nt);
if isa(u.value,'MULTIMATRIX')
   v{k} = reshape(u.value{k},u.s);
elseif isa(u.value,'double')
   for k=1:nt
    v{k} = reshape(u.value(:,k),u.s);
   end
else
error('pas prevu')
end
u.value = v;

end