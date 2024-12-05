function w = sumelem(u)
% function w = sumelem(u)

if ~israndom(u)
    for k=1:length(u.value)
        v{k} = sum(u.value{k},3);
        v{k} = double(v{k});
    end
    w = v{1};
    for k=2:length(v)
        w = w+v{k};
    end
else
    for k=1:length(u.value)
        v{k} = sum(u.value{k},3);
        v{k} = double(v{k});
    end
    w = v{1};
    for k=2:length(v)
        w = w+v{k};
    end
end
