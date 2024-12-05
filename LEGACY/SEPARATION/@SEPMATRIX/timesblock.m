function M = timesblock(A)

i=1;
Mi = A.F{i,1};
for k=2:A.dim
Mi = Mi.*A.F{i,k};    
end
M = Mi*A.alpha(i);

for i=2:A.m
Mi = A.F{i,1}; 
for k=2:A.dim
Mi = Mi.*A.F{i,k};    
end
M = M + Mi*A.alpha(i);
end

return