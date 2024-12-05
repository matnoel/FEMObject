function y = SEPVECTOR(x)

y = SEPVECTOR(x.funs{1});
for i=2:length(x.funs)
   y = y + SEPVECTOR(x.funs{i}); 
end
