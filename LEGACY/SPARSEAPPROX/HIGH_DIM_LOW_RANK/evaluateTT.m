function fs=evaluateTT(v,Hs)
n= size(Hs{1},1);
d=numel(Hs);
fsTest=cell(n,1);
fsTest(:)={1}; 
A=cell(1,d);
    for j=1:d
        A{j}=ttm(v{j},Hs{j},3);
    end
fs=ones(n,1);
for mu=1:d
    fs = repmat(fs,[1 1 size(A{mu},2)]);
    fs = sum(fs.*permute(A{mu},[3 1 2]),2);
    fs = fs(:,:);
end

