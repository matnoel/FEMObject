function points = calc_collocationpoints(PC)

points = calc_collocationpoints(PC.RANDPOLYS,PC.indices(:,1:end-1)+1);
