function stodim = getstodim(rv)

for i=1:length(rv.Y)
    stodim{i}=getnumber(rv.Y{i});
end

stodim = [stodim,{rv.number}];