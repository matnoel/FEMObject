function A = std(u)

for k=1:u.M
A{k}=std(u.RV{k});
end
