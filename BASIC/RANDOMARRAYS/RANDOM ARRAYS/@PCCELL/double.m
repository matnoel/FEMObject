function apc=double(apc)

for i=0:getP(apc)
   apc.value{i+1}=double(apc.value{i+1}); 
end