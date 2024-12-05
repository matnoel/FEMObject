function ok = eq(r1,r2)
% function ok = eq(r1,r2)

if ~(isa(r1,'RANDVARS') && isa(r2,'RANDVARS') )
    ok=0;
elseif r1.M~=r2.M
    ok=0;
else
    ok=1;
    for i=1:r1.M
        ok = ok & (r1.RV{i} == r2.RV{i});
    end
end