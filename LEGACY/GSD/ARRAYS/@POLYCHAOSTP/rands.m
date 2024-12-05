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

Hm = cell(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
    Hm{i} = rand(length(PC,i),1);
end

x = PCTPMATRIX(PC,rand(s),Hm{:});