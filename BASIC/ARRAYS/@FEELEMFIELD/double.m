function v = double(u)
% function v = double(u)

if length(u.value)==1
    v = double(u.value{1});
    for i=2:length(u.value)
        v = cat(3,v,double(u.value{i}));
    end
else
    error('double ne peut pas marcher : plusieurs groupes d''elements')
end