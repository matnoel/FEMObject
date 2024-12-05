function size = subplot_format(n)
% size = subplot_format(n)
% Determine suitable subplot array format according to number of plots
c = ceil(sqrt(n)) ;
l = ceil(n/c) ;
size = [l c] ;
end
