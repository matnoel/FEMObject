function x = calc_midpoint(elem,node)

x = calc_midpointelem(elem,node);

x=POINT(mean(double(x),1));