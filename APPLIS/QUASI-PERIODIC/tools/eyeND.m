function I = eyeND(sz)
subs = repmat((1:min(sz))',1,numel(sz)) ;
I = zeros(sz) ;
I(ind) = 1 ;
end