function gauss = elem_gauss(elem,ordre)
% function gauss = elem_gauss(elem,ordre)

npoints = ceil((ordre+1)/2);
gauss = calc_gausspoints(POLYLEGENDRE(),npoints);
gauss.w = gauss.w*2;
