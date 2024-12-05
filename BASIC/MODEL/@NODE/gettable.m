function table = gettable(u)

if ~isempty(u.number)
table=[u.number,double(getcoord(u))];
else
table=zeros(0,size(double(getcoord(u)),2)+1);    
end