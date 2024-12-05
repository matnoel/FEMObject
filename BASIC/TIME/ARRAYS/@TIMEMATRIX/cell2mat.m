function u=cell2mat(u)

if iscell(u)
for k=1:length(u.value)
    if isa(u.value{k},'double')
   u.value{k} = reshape(u.value{k},prod(u.s),1) ;
    else
        error('pas le bon type')
    end
end
    u.value = [u.value{:}];
end