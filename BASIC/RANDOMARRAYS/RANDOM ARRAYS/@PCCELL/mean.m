function am=mean(apc)

Mfe = findpolyfe(apc);

if isempty(Mfe)
am = apc.value{1};

else

Hm = mean(getPC(apc));
rep = find(Hm~=0);
am = apc.value{rep(1)}*Hm(rep(1));
for k=2:length(rep) 
 am = am + apc.value{rep(k)}*Hm(rep(k));  
end

end