function gauss = elem_gauss(elem,order)
% function gauss = elem_gauss(elem,order)

npoints = ceil((order+1)/2);
gauss = calc_gausspoints(RANDPOLYS(POLYLEGENDRE(),POLYLEGENDRE(),POLYLEGENDRE()),npoints);
gauss.w = gauss.w*8;
