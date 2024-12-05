function PC = calc_mean(PC)

Hm = cell(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
  Hm{i} = mean(PC.PCgroups{i});
end

PC.mean = Hm;

