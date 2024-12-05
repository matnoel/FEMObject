function metricuni = calc_metricuni(h,p,varargin)

metricuni=cell(1,h.M);
for k=1:h.M
metricuni{k}=calc_metric(h.h{k},p(k),varargin{:});
end
