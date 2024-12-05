function Z=sparse(Z)

for i=1:length(Z.value)
  Z.value{i}=sparse(Z.value{i});  
end