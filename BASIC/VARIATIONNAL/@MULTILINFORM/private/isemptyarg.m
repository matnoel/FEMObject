function r = isemptyarg(arg)
r=[];
for i=1:length(arg)
r=[r, isempty(arg{i})];
end

