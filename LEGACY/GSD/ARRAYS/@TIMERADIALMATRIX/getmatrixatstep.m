function Z = getmatrixatstep(atr,i)

if isa(atr.V,'FEELEMFIELD')
Z = atr.V * atr.D * double(getmatrixatstep(atr.L,i))  ;
elseif iscell(atr.V)
Z = multimtimes(atr.V,atr.D*double(getmatrixatstep(atr.L,i)));    
Z = double(cell2mat(Z));
Z = reshape(Z,size(atr));
else
Z = double(atr.V)* atr.D * double(getmatrixatstep(atr.L,i)) ;   
Z = reshape(Z,size(atr));
end
