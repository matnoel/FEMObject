function m = get_ximasse(x,i)

if isempty(x.ximasse)
    error('ximasse non calculee')
end

if nargin==1
    m = x.ximasse;
elseif length(i)==1
    m = x.ximasse{i};
else
    m = x.ximasse(i);
end

