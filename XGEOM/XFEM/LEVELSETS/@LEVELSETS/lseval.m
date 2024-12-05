function ls=lseval(ls,P,varargin)
 
for i=1:ls.n
   LS{i} = lseval(ls.LS{i},P,varargin{:}) ;
end

if isa(P,'MODEL') | ischarin('LEVELSET',varargin)
   ls.LS = LS ;
else
   ls = LS;
end
