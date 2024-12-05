function PC = calc_metricgroups(PC,varargin)

for k=1:length(PC.groups)
   PC.PCgroups{k} = calc_metric(PC.PCgroups{k});
end