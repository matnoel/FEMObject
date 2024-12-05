function x = transpose(x)

for i=1:length(x.funs)
x.funs{i} = transpose(x.funs{i});
end
