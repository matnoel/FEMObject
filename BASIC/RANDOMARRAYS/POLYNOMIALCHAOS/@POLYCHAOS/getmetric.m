function masse=getmetric(PC)

masse = PC.metric;

if isempty(masse)
    error('metric not computed : use calc_metric')
end

