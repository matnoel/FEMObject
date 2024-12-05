function n=getnbelem(H)

n = 1;
for k=1:getM(H)
n=n * getnbelem(H.h{k});   
end

