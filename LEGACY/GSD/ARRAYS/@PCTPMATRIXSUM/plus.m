function x = plus(x,y)

if isa(x,'PCTPMATRIXSUM') && isa(y,'PCTPMATRIXSUM')
    x.funs = [x.funs,y.funs];
elseif isa(x,'PCTPMATRIXSUM') && isa(y,'PCTPMATRIX') 
    x.funs = [x.funs,{y}];
elseif isa(y,'PCTPMATRIXSUM') && isa(x,'PCTPMATRIX')    
    y.funs = [{x},y.funs];
    x=y;
elseif isa(x,'PCTPMATRIXSUM') && isa(y,'double')
    y = PCTPMATRIX(x.POLYCHAOSTP,y);
    x = plus(x,y);
elseif isa(y,'PCTPMATRIXSUM') && isa(x,'double')   
    x = PCTPMATRIX(y.POLYCHAOSTP,x);
    x = plus(x,y);
else
    error('pas prevu')
end