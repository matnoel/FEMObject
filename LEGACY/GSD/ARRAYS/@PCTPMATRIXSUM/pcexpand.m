function apc = pcexpand(a)

if ~isempty(a.funs)
apc = pcexpand(a.funs{1});
for i=2:length(a.funs)
    apc = apc + pcexpand(a.funs{i});
end
else
   warning('aucune fonction : pcexpand non realise') 
end

    