function u=reshapem(u,a,b)
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

    