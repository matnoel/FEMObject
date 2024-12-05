function Z=expand(atr)

if iscell(atr.V)
    atr.V = cell2mat(atr.V);
    Z = expand(atr);
else
    Z = double(atr.V)* atr.D * double(atr.L) ;
    Z = TIMEMATRIX(Z,atr.TIMEMODEL,size(atr));
end

