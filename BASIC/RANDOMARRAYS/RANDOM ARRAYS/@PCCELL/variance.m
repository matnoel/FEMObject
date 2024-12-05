function av=variance(apc)

av=apc.value{1};
for k=2:apc.value 
 av = av + apc.value{k}.^2;  
end
av=av-mean(apc).^2;
