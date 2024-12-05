function u=reshapem(u,a,b)

if isa(u.value,'cell')
if nargin==2
    u.value = reshape(u.value,a);
else
    u.value = reshape(u.value,a,b);    
end
u.sm = size(u.value);

else
    
if nargin==3
  if ~(a*b == prod(u.sm))
      error('To RESHAPE the number of elements must not change.')
  end
  u.sm=[a,b];
elseif nargin==2

  if ~(prod(a) == prod(u.sm))
      error('To RESHAPE the number of elements must not change.')
  end
  u.sm = a;
else
    error('enter 2 or 3 arguments')
end

    
end