function a = cellctranspose(a)
%function a = cellctranspose(a)

a = cellfun(@ctranspose,a,'UniformOutput',false);
%for k=1:numel(a)
%   a{k} = ctranspose(a{k}); 
%end
