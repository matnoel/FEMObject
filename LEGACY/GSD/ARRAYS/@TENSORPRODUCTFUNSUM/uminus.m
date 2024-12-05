function a = uminus(a)
% function a = uminus(a)

for i=1:length(a.tensorfuns) 
    a.tensorfuns{i}=uminus(a.tensorfuns{i});
end
