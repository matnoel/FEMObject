function rpc = rands(varargin)
% function rpc = rands(varargin)

PC=POLYCHAOS(getclassin('POLYCHAOS',varargin));
s = getclassin('double',varargin);
if isa(s,'cell')
    s = [s{:}];
end

if length(s)==1
    s=[s,s];
end

rpc = rand(prod(s),length(PC));
rpc = PCMATRIX(rpc,s,PC);
