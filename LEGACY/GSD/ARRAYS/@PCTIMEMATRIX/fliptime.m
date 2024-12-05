function u = fliptime(u)

if isa(u.value,'cell')
u.value = u.value(end:-1:1);

elseif isa(u.value,'PCMATRIX') 

u.value = u.value(:,length(u.TIMEMODEL):-1:1);    
    
elseif isa(u.value,'PCRADIALMATRIX')
    
u.value = u.value(:,length(u.TIMEMODEL):-1:1);    
 
else
    error('pas prevu');
end
