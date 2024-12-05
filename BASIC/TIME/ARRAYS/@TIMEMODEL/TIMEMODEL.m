function T = TIMEMODEL(t0,t1,n)
% function T = TIMEMODEL(t0,t1,n)
% function T = TIMEMODEL(t)

switch nargin
case 0
    t = [0,1];
    t0 = 0;
    t1 = 1;
case 1
    t = t0;
    t0 = t(1);
    t1 = t(end);
case 3
    t = linspace(t0,t1,n+1);
otherwise
    error('rentrer un ou 3 parametres')
end

dt = t(2:end)-t(1:end-1);
nt = length(dt);

T.t = t;
T.t0 = t0;
T.t1 = t1;
T.dt = dt;
T.nt = nt;
T.uniform = isuniform(T.dt);
T.approxparam = PARAMETERS('p',1,'type','default');
T.evolparam = PLOTOPTIONS();
T.evolparam = setparam(T.evolparam,'makepause',true,'pausetime',1/nt,'plotiter',false,'plottime',true,'plotstep',1);

T = class(T,'TIMEMODEL');
