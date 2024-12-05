function masseuni = calc_masseuni(h,p,p2,varargin)
if nargin==2 || isempty(p2)
p2=p;
end

masseuni=cell(1,h.M);
for k=1:h.M
masseuni{k}=calc_masse(h.h{k},p(k),p2(k),varargin{:});
end
