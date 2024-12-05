function a = squeeze(a)
% function x = PCMATRIX(a,s,PC)

n = length(a.s);
rep=[];
for k=n:-1:3
    if a.s(k)==1
        rep=[rep,k];
    else
        break
    end
end
a.s(rep)=[];