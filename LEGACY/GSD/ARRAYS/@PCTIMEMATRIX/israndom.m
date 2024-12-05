function rep = israndom(u)
% function rep = israndom(u)

if isa(u.value,'cell')
    rep = false;
    for k=1:length(u.value)
        rep = rep | israndom(u.value{k});
        if rep
            break
        end
    end
    
else
    rep = israndom(u.value);
end
