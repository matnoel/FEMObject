function gauss=calc_gausspoints(h,n)
% function gauss=calc_gausspoints(h,n)
% calcul des points de gauss par morceau pour une base FE sur [0,1]
%

gauss = calc_gausspoints(h.POLYFE,n);
gauss.coord = transfer(RANDVAR(h.POLYFE),h.rv,gauss.coord);
