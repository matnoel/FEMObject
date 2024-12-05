function B = PCTPMATRIX(A,PC,reppc,varargin)
% function B = PCTPMATRIX(A,PC,reppc)

noexpand=getcharin('noexpand',varargin,0);

dim = getdim(A);
dimpc = getnbgroups(PC);
if nargin==2
    reppc = [dim-dimpc+1:dim];
end
repnopc = setdiff(1:getdim(A),reppc);

if isempty(repnopc)
    i=1;
    B = PCTPMATRIX(PC,A.alpha(i),A.F{i,reppc});
    for i=2:getm(A)
        B = B + PCTPMATRIX(PC,A.alpha(i),A.F{i,reppc});
    end
else
    tempA=removedim(A,reppc);
    if noexpand
        temp=truncate(tempA,1);
    else
        temp=expand(truncate(tempA,1));
    end
    B = PCTPMATRIX(PC,temp,A.F{1,reppc});
    for i=2:getm(A)
         if noexpand
            temp=truncate(tempA,i);
         else
            temp=expand(truncate(tempA,i));
         end
        B = B + PCTPMATRIX(PC,temp,A.F{i,reppc});
    end
end


