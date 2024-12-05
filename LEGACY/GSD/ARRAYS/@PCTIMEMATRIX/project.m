function u = project(u,PC)
% function u = project(u,PC)

if isa(u.value,'cell')
    for k=1:length(u.value)
        u.value{k}=project(u.value{k},PC);
    end
else
    u.value=project(u.value,PC);
end
