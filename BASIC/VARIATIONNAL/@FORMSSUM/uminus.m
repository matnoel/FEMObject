function a = uminus(a)
% function a = uminus(a)

for k=1:length(a.forms)
    a.forms{k} = uminus(a.forms{k});
end