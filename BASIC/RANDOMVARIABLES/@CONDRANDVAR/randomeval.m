function X = randomeval(rv,y,RV)
% function X = randomeval(rv,y,RV)

alleval = 1;
for i=1:length(rv.Y)
    if isa(rv.Y{i},'RANDVAR') || isa(rv.Y{i},'CONDRANDVAR')
        if ismember(rv.Y{i},RV)
            rv.Y{i} = randomeval(rv.Y{i},y,RV);
        else
            alleval = 0;
        end
    end
end

if alleval
    for i=1:length(rv.funparam)
        if isa(rv.funparam{i},'inline') || isa(rv.funparam{i},'function_handle')
            rv.funparam{i} = rv.funparam{i}(rv.Y{:});
        end
    end
    X = setnumber(rv.X(rv.funparam{:}),rv.number);
    
    [ok,rep] = ismember(X,RV);
    if ok
        X = randomeval(X,y,RV);
    end
else
    X = rv ;
end
