function a = celltranspose(a)
%function a = celltranspose(a)

for k=1:numel(a)
    a{k} = transpose(a{k}); 
end
