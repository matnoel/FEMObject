function onepc=ones(varargin)
% function onepc=ones(varargin)

PC = POLYCHAOS(getclassin('POLYCHAOS',varargin));
s = getclassin('double',varargin);
if isa(s,'cell')
    s = [s{:}];
end

if length(s)==1
    s = [s,s];
end

onepc = full(mean(PC)');
onepc = repmat(onepc,prod(s),1);
onepc = PCMATRIX(onepc,s,PC);
