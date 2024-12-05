function gauss = elem_gauss(elem,order)
% function gauss = elem_gauss(elem,order)

npoints = ceil((order+1)/2);
gauss = calc_gausspoints(POLYLEGENDRE(),npoints);
gauss.w = gauss.w*2;
