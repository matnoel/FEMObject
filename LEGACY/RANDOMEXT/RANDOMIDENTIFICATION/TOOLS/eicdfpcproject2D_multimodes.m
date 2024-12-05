function [Xpc,Ym,P] = eicdfpcproject2D_multimodes(Xblock,p,h1,h2,varargin)
% function [Xpc,Ym,P] = eicdfpcproject(Xblock,P,p,h1,h2,z)
% p : ordre du chaos
% h1,h2: RANDPOLY
% 
% See also eicdfpcproject, separatesamples

m=length(Xblock);
Q = 0;
for l=1:m
    Q = Q+size(Xblock{l},2);
end
P = zeros(1,m);
x = zeros(1,m+1);
for l=1:length(Xblock)
    P(l) = size(Xblock{l},2)/Q;
    x(l+1) = x(l)+P(l);
end

if nargin<=2 || isempty(h1)
    h1 = POLYHERMITE();
end
if nargin<=3 || isempty(h2)
    h2 = POLYHERMITE();
end

PC = POLYCHAOS(RANDPOLYS(POLYFE(x),h1,h2),[0,p,p],'typebase',2);
PC = restrictorder(PC,p,2:3,1);

ind = getindices(PC); 
Xpc = zeros(2,length(PC));
Ym=cell(1,m);

for i=1:m
    Ym{i} = eicdfpcproject2D(Xblock{i},p,h1,h2,varargin{:});
    rep = find(ind(:,1)==i-1);
    Xpc(:,rep)=sqrt(P(i))*double(Ym{i});
end

Xpc = PCMATRIX(Xpc,[2,1],PC);






