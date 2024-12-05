function myparallel(cmd,varargin)
% function myparallel(cmd,varargin)
% Parallel Computing Toolbox
% cmd : char ('start' or 'open', 'stop' or 'close', 'find' or 'delete')
%
% myparallel('start')
% myparallel('open')
% start a parallel pool with one worker per core available on the local machine, up to 12
%
% myparallel('start',nbWorkers)
% myparallel('open',nbWorkers)
% start a parallel pool with nbWorkers local workers requested (up to 4 in R2006b, increased to 8 in R2009a, and to 12 in R2011b); 
% it should not exceed the number of cores of the local machine c.NumWorkers
%
% myparallel('stop')
% myparallel('close')
% shut down and delete the parallel pool
%
% myparallel('find')
% find all job objects stored in cluster and sorted by state
%
% myparallel('find',state)
% state : char ('pending', 'queued', 'running' or 'finished')
% find all job or task objects whose State is state fom cluster
%
% myparallel('find',i)
% i : double
% find the job object identified by Id i on cluster
%
% myparallel('delete')
% delete all job or task objects from cluster
%
% myparallel('delete',state)
% state : char ('pending', 'queued', 'running' or 'finished')
% delete all job or task objects whose State is state fom cluster
%
% myparallel('delete',i)
% i : double
% delete the job object identified by Id i from cluster
%

if ~usejava('jvm')
    warning('The Java Virtual Machine is not supported in MATLAB.');
    return
end
if isempty(ver('parallel'))
    warning('The Parallel Computing Toolbox is not installed.');
    return
end
if ~license('test','Distrib_Computing_Toolbox')
    warning('No license exists for the Parallel Computing Toolbox.');
    return
end
if ~license('checkout','Distrib_Computing_Toolbox')
    warning('No license checked out for the Parallel Computing Toolbox.');
    return
end

p = gcp('nocreate'); % returns the current pool if one exists
c = parcluster; % returns a cluster object representing the cluster identified by the default cluster profile
if nargin==0
    if isempty(p)
        cmd = 'start';
    else
        cmd = 'stop';
    end
end
switch cmd
    case 'disp'
        disp(c)
    case {'start','open'}
        if isempty(p)
            if nargin<2
                % matlabpool(c); % open a parallel pool with one worker per core available on the local machine, up to 12
                parpool(c); % new command from R2013b
            else
                nbWorkers = varargin{1}; % set the number of local workers
                % matlabpool(c,nbWorkers); % open a parallel pool with nbWorkers local workers
                parpool(c,nbWorkers); % new command from R2013b
            end
        end
        % fprintf(['\nNumber of currently open workers  = ' num2str(PoolSize) '\r']);
    case {'stop','close'}
        if ~isempty(p)
            % matlabpool close % close the parallel pool
            delete(p) % new command from R2013b
        end
    case 'find'
        if nargin<2
            [pending queued running finished] = findJob(c)
        else
            if ischar(varargin{1})
                findJob(c,'State',varargin{1})
            else
                findJob(c,'Id',varargin{1})
            end
        end
    case 'delete'
        if nargin<2
            delete(c.Jobs)
        else
            if ischar(varargin{1})
                j = findJob(c,'State',varargin{1});
            else
                j = findJob(c,'Id',varargin{1});
            end
            delete(j)
        end
    otherwise
        error('Wrong command');
end

end
