function B = PCTPRADIALMATRIX(A,PC,reppc)
dim = getdim(A);
dimpc = getnbgroups(PC);
if nargin==2
reppc = [dim-dimpc+1:dim];
end
repnopc = setdiff(1:getdim(A),reppc);

if length(repnopc)==1
V = extractvectors(A,repnopc)*diag(A.alpha);
i=1;
L = PCTPMATRIX(PC,getei(i,A.m),A.F{i,reppc});
for i=1:A.m
L = L + PCTPMATRIX(PC,getei(i,A.m),A.F{i,reppc});
end

B = PCTPRADIALMATRIX()

else
    error('pas programme')
end


function u = getei(i,n)

u = zeros(n,1);
u(i)=1;

