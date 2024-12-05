function x = uminus(x)
% function x = uminus(x)

for i=1:length(x.funs)
   x.funs{i} = uminus(x.funs{i}); 
end