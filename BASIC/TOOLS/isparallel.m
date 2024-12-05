function [isPoolOpen,PoolSize] = isparallel
% function [isPoolOpen,PoolSize] = isparallel
% Parallel Computing Toolbox
% Check if a parallel pool is open and return the size of the parallel pool

if ~usejava('jvm')
    % The Java Virtual Machine is not supported in MATLAB
    PoolSize = 0;
    isPoolOpen = false;
else
    PoolSize = matlabpool('size');
    isPoolOpen = (PoolSize > 0);
end

end
