function n = normbyelement(u)

PC = getPC(u);
n = zeros(1,getnbelem(PC));
v = double(u);
for i=1:getnbelem(PC)
rep = getindexofelement(PC,i);
t = v(:,rep);
n(i) = sqrt(sum(sum(t.^2)));
end


