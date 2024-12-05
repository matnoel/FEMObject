function ok = eq(r1,r2)
% function ok = eq(r1,r2)

if ~(isa(r1,'CONDRANDVAR') && isa(r2,'CONDRANDVAR') )
    ok = 0;
elseif ~strcmp(func2str(r1.X),func2str(r2.X))
    ok = 0;
elseif length(r1.Y)~=length(r2.Y)
    ok = 0;
else
    ok = 1;
    for i=1:length(r1.Y)
        ok = ok & (r1.Y{i} == r2.Y{i});
    end
end