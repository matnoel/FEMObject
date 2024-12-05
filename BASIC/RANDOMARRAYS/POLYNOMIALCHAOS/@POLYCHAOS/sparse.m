function zeropc=sparse(varargin)

PC = getclassin('POLYCHAOS',varargin);
s=getclassin('double',varargin);

if isa(s,'cell')
    s=[s{:}];
end
if length(s)==1
    s=[s,1];
end

zeropc = PCMATRIX(sparse(prod(s),length(PC)),s,PC);