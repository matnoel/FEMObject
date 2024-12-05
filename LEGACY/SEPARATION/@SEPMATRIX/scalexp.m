function A=scalexp(A)
% H=scalexp(H)
% Certaines manipulations (evaluation/random) reduisent des dimension
% a un scalaire : il faut faire un peu de menage

s=size(A);
scaldim=find(s==1);
vectdim=find(s~=1);
if ~isempty(scaldim)
    A.alpha=prod([prod(full(cell2mat(A.F(:,scaldim))),2)';A.alpha],1);
else
    return
end
if ~isempty(vectdim)
    A.F=A.F(:,vectdim);
    A.dim=A.dim-length(scaldim);
else
    A=sum(A.alpha);
end
