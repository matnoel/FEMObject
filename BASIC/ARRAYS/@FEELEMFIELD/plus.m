function s = plus(s1,s2)
% function s = plus(s1,s2)

if isa(s1,'FEELEMFIELD') && isa(s2,'FEELEMFIELD')
    try
        s = s1;
        for k=1:length(s1.value)
            s.value{k} = s1.value{k}+s2.value{k};
        end
    catch
        error('les FEELEMFIELD ne concordent pas')
    end
elseif isa(s1,'FEELEMFIELD')
    s = s1;
    for k=1:length(s1.value)
        s.value{k} = s1.value{k}+s2;
    end
elseif isa(s2,'FEELEMFIELD')
    s = s2;
    for k=1:length(s2.value)
        s.value{k} = s2.value{k}+s1;
    end
end
