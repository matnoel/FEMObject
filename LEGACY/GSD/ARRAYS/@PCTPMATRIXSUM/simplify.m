function x = simplify(x)

isr = isranddim(x);
rdim = find(isr);
if isempty(rdim) 
y = simplify(x.funs{1});
for i=2:length(x.funs)
y = y +  simplify(x.funs{i});   
end
x=y;
elseif length(rdim)==1
y = 0;
for i=1:length(x.funs)   
if isranddim(x.funs{i},rdim)    
y = y + simplify(x.funs{i});
else
y = y + mean(x.POLYCHAOSTP,rdim)*simplify(x.funs{i});
end
end
x=y;
end
