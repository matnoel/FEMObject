function y = SEPMATRIX(x)

y = SEPMATRIX(x.funs{1});
for i=2:length(x.funs)
   y = y + SEPMATRIX(x.funs{i}); 
end
