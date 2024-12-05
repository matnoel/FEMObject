function u = PCMATRIX(u)
T = u.TIMEMODEL;
nt = length(T);
if isa(u.value,'cell')
PC = getPC(u);
s = u.s ; 
v = cell(1,numel(PC));

for k=1:numel(PC)
v{k} = getpccompo(u,k);
v{k} = getvalue(cell2mat(v{k}));
end

u = PCMATRIX(v,[prod(u.s),nt],PC);

elseif isa(u.value,'PCMATRIX') |  isa(u.value,'PCRADIALMATRIX')
    
u = reshape(u.value,prod(u.s),nt); 
    
else
    error('conversion non prevue');
end
