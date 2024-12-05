function [rs,fs,RV] = orange(numSample)
RV = RANDVARS(RVUNIFORM(-1,1),4);
n = random(RV,1,numSample);
rs = (cell2mat(n'))';
fs=zeros(numSample,1);
for i=1:numSample
    fs(i,1) = metamodele(rs(i,:));
end
end