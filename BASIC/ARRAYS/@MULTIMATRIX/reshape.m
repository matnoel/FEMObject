function u=reshape(u,a,b)

if isa(u.value,'cell')
    if nargin==2
   for k=1:numel(u.value)
      u.value{k}=reshape(u.value{k},a); 
   end
   u.s = a ; 
    elseif nargin==3
        for k=1:numel(u.value)
      u.value{k}=reshape(u.value{k},a,b); 
        end 
   u.s = [a,b];
    end

else
if nargin==3
  if ~(a*b == prod(u.s))
      error('To RESHAPE the number of elements must not change.')
  end
  u.s=[a,b];
elseif nargin==2
    if ~(prod(a) == prod(u.s))
      error('To RESHAPE the number of elements must not change.')
  end
  u.s = a;
else
    error('enter 2 or 3 arguments')
end


end