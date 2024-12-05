function zeropc = zeros(varargin)
% function zeropc = zeros(varargin)

PC = getclassin('POLYCHAOS',varargin);
s=getclassin('double',varargin);

if isa(s,'cell')
    s=[s{:}];
end
if length(s)==1
    s=[s,1];
end

zeropc = PCMATRIX(zeros(prod(s),length(PC)),s,PC);