function u = reshape(u,s1,s2)

if nargin==2 && length(s1)>=2
    if prod(s1)~=prod(u.s)
      error('To reshape the number of elements must not change')
    end
    u.s = s1;
elseif nargin==3 && length(s1)==1 && length(s2)==1
     if (s1*s2)~=prod(u.s)
      error('To reshape the number of elements must not change')
     end
    u.s=[s1,s2];
else
    error('mauvais arguments')
end