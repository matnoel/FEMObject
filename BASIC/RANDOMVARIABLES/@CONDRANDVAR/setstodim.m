function rv = setstodim(rv,stodim)
if isa(stodim,'cell')
stodim = [stodim{:}];
end

if length(stodim)~=length(rv.Y)+1
    error('rentrer la dimension stochastique pour toutes les variables')
end

for i=1:length(rv.Y)
if isa(rv.Y{i},'RANDVAR')
    rv.Y{i}=setnumber(rv.Y{i},stodim(i));
elseif isa(rv.Y{i},'CONDRANDVAR')
    rv.Y{i}=setstodim(rv.Y{i},stodim(1:i));
end
end

rv.number = stodim(end);

