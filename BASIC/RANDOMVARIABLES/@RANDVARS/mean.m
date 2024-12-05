function A = mean(u)

for k=1:u.M
A{k}=mean(u.RV{k});
end
