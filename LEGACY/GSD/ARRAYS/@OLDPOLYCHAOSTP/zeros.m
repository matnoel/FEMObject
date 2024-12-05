function x = zeros(varargin)
% function x = zeros(varargin)

PC = getclassin('POLYCHAOSTP',varargin);
s=getclassin('double',varargin);

if isa(s,'cell')
    s=[s{:}];
end
if length(s)==1
    s=[s,1];
end

Hm = cell(1,getnbdim(PC));
for i=1:getnbdim(PC)
    Hm{i} = zeros(getn(PC,i),1);
end

x = PCTPMATRIX(PC,zeros(s),Hm{:});