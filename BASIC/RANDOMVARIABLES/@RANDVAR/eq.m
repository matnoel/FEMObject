function ok = eq(r1,r2)
% function ok = eq(r1,r2)

if ~(isa(r1,'RANDVAR') && isa(r2,'RANDVAR') )
    ok = 0;
elseif ~strcmp(class(r1),class(r2))
    ok = 0;
else
    ok = 1;
    for i=1:size(r1.param,1)
        ok = ok & (r1.param{i,2} == r2.param{i,2});
    end
end