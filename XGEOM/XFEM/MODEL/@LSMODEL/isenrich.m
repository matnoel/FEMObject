function rep = isenrich(S)

rep=false;
for i=1:getnbgroupelem(S)
   if isenrich(getgroupelem(S,i))
       rep=true;
       break
   end
end
