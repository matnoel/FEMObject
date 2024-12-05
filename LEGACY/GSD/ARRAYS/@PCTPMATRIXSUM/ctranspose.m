function x = ctranspose(x)

for i=1:length(x.funs)
x.funs{i} = ctranspose(x.funs{i});
end
