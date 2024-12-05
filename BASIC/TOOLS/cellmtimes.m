function c = cellmtimes(a,b)%,funa,funb)
%function c = cellmtimes(a,b)


if ~isa(a,'cell') && ~isa(b,'cell')
    error('un argument doit etre une cellule')
elseif isa(a,'cell') && ~isa(b,'cell')
    c=a;
    for k=1:numel(a)
        c{k} = a{k}*b;
    end    
elseif isa(b,'cell') && ~isa(a,'cell')
    c=b;
    for k=1:numel(b)
        c{k} = a*b{k};
    end   
elseif isa(b,'cell') && isa(a,'cell') 
    if ~all(size(b)==size(a))
        error('le nombre de cellules ne correspond pas')
    end
    c= a ;
    for k=1:numel(a)
        c{k} = a{k}*b{k};    
    end
end




