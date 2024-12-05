function PC = calc_mean(PC)

Hm = cell(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
  Hm{i} = mean(getpoly(PC,i),0:getn(PC,i)-1);
end

PC.mean = Hm;

