function PC=calc_masseuni(PC,PC2,varargin)
p=PC.p;
if nargin==1 || isempty(PC2)
    p2 = p;
elseif isa(PC2,'POLYCHAOS')
    p2 = getorder(PC2);
end

PC.masseuni = calc_masseuni(PC.RANDPOLYS,p,p2,varargin{:});
