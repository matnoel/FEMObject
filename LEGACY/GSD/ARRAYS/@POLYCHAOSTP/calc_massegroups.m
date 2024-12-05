function PC = calc_massegroups(PC,varargin)

for k=1:length(PC.groups)
   PC.PCgroups{k} = calc_masse(PC.PCgroups{k});
end