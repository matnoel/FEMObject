function [fsReg,AReg,fsTest,ATest]=cross_validation_Kfold_least_squares(A,fs,K)
fsReg=cell(K,1);AReg=cell(K,1); fsTest=cell(K,1); ATest=cell(K,1);
CVO=cvpartition(length(fs),'kfold',K);
for i=1:K
    trIdx=CVO.training(i);
    if ndims(A)==2
        fsReg{i}=fs(trIdx);AReg{i}=A(trIdx,:);
        fsTest{i}=fs(~trIdx); ATest{i}=A(~trIdx,:);
    elseif ndims(A)==3
        fsReg{i}=fs(trIdx);AReg{i}=A(:,:,trIdx);
        fsTest{i}=fs(~trIdx); ATest{i}=A(:,:,~trIdx,:);  
    end
end