function u=reshape(u,a,b)
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

