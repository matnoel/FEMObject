function PC=calc_metricuni(PC,varargin)
p=PC.p ;

PC.metricuni = calc_metricuni(PC.RANDPOLYS,p,varargin{:});
