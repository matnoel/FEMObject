function param = randomsetparam(param,RVvalue,RV)
% function param = randomsetparam(param,RVvalue,RV)

for i=1:size(param,1)
    if isa(param{i,2},'RANDVAR')
        [ok,rep] = ismember(param{i,2},RV);
        if ok
            param{i,2} = RVvalue{rep};
        end
    end
    
end
