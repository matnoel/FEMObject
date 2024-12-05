function apc = project(a,pc,varargin)
% function apc = project(a,pc,varargin)

apc = cell(1,a.M);
for i=1:a.M
    apc{i} = project(a.RV{i},pc,varargin{:});
end
apc = vertcat(apc{:});
if a.M==0
    apc = PCMATRIX(zeros(0,1),[1,0],pc);
end
