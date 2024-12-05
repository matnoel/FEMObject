function x = mtimes(x,y)

if isa(y,'PCTPMATRIXSUM') && ~isa(x,'PCTPMATRIXSUM')    
    
    for i=1:length(y.funs)
        y.funs{i} = mtimes(x,y.funs{i});
    end
    x=y;
elseif isa(x,'PCTPMATRIXSUM') && ~isa(y,'PCTPMATRIXSUM') 
  
    for i=1:length(x.funs)
        x.funs{i} = mtimes(x.funs{i},y);
    end
elseif isa(x,'PCTPMATRIXSUM') && isa(y,'PCTPMATRIXSUM') 
    z = PCTPMATRIXSUM(x.POLYCHAOSTP);
    z.funs = cell(1,length(x.funs)*length(y.funs));
    ll=0;
    for i=1:length(x.funs)
        for j=1:length(y.funs)
            ll=ll+1;
        z.funs{ll} = mtimes(x.funs{i},y.funs{j});
        end
    end
else
    error('pas prevu')
end