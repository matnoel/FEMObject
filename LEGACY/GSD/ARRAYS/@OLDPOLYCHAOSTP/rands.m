function x = rands(varargin)
% function x = rands(varargin)

PC=POLYCHAOSTP(getclassin('POLYCHAOSTP',varargin));
s = getclassin('double',varargin);
if isa(s,'cell')
    s = [s{:}];
end

if length(s)==1
    s=[s,s];
end

Hm = cell(1,getnbdim(PC));
for i=1:getnbdim(PC)
    Hm{i} = rand(getn(PC,i),1);
end

x = PCTPMATRIX(PC,rand(s),Hm{:});