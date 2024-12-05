function rv = sort(rv)

number=getnumber(rv);
number=[number{:}];
[a,b]=sort(number);
rv.RV = rv.RV(b);
