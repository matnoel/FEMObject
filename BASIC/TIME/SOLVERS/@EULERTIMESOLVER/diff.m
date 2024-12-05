function vt = diff(L,ut)
% function vt = diff(L,ut)
% Computes v = u'
% L: EULERTIMESOLVER
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'

eulertype = getparam(L,'eulertype');

T = L.TIMEMODEL;
% t = gett(T);
nt = getnt(T);
dt = getdt(T);

n = size(ut,1);
vt = cell(1,length(T));

switch eulertype
    case 'explicit'
        for i=1:nt
            vt{i} = (getmatrixatstep(ut,i+1) - getmatrixatstep(ut,i))/dt(i);
        end
        vt{nt+1} = zeros(n,1);
    case 'implicit'
        vt{1} = zeros(n,1);
        for i=1:nt
            vt{i+1} = (getmatrixatstep(ut,i+1) - getmatrixatstep(ut,i))/dt(i);
        end
end

vt = TIMEMATRIX(vt,T,[size(ut,1),1]);
