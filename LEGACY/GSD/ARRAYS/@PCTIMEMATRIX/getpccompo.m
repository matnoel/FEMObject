function v = getpccompo(u,k)
T = u.TIMEMODEL;
if isa(u.value,'cell')
v = cell(1,length(u.value));
for i=1:length(u.value)
v{i} = getpccompo(u.value{i},k);
end
v = TIMEMATRIX(v,T,u.s);

elseif isa(u.value,'PCMATRIX') || isa(u.value,'PCRADIALMATRIX')
v = getpccompo(u.value,k);
if all(u.s>1)
    error('pas prevu')
end
v = reshape(v,prod(u.s),length(gettapprox(T)));
v = TIMEMATRIX(v,T,u.s);  
 
else
    error('conversion non prevue');
end
