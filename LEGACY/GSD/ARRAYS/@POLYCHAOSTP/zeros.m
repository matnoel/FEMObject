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
for i=1:getnbgroups(PC)
    Hm{i} = zeros(length(PC.PCgroups{i}),1);
end

x = PCTPMATRIX(PC,ones(s),Hm{:});