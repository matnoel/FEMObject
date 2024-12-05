function onepc=ones(varargin)

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
  Hm{i} = mean(PC,i);%full(double(mean(getpoly(PC,i),0:getn(PC,i)-1)));   
end

onepc = PCTPMATRIX(PC,ones(s),Hm{:});