function rep=isuniform(dt)
% function isuniform(dt)
% teste si les pas de temps sont constants

dt1 = dt(1);

if any(abs(dt-dt1)>10*eps)
    rep=false;
else
    rep=true;
end
