function x = plus(x,y)

if isa(x,'PCTPMATRIX') && isa(y,'PCTPMATRIX')
    
    x = PCTPMATRIXSUM(x,y);
    
elseif isa(x,'double')
    
    x = PCTPMATRIXSUM(x*one(y.POLYCHAOSTP),y);
    
elseif isa(y,'double')
    
   x = PCTPMATRIXSUM(x,y*one(x.POLYCHAOSTP));
    
else
    error('pas programme')
    
end