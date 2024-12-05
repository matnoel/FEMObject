function [U,w,lambda] = initializeVariables(glob,patches,interfaces)
% function [U,w,lambda] = initializeVariables(glob,patches,interfaces)
% Initializes global solution U, local solution w and Lagrange multipier lambda
%
% Inputs:
% glob: Global
% patches: Patches
% interfaces: Interfaces
%
% Outputs:
% U: initial global solution U
% w: initial local solution w
% lambda: initial Lagrange multiplier lambda

if isempty(glob.timeSolver)
    U = zeros(getnbddlfree(glob.S),1);
else
    T = gettimemodel(glob.timeSolver);
    sz = [getnbddlfree(glob.S),getnbtimedof(glob.timeSolver)];
    U = zeros(sz);
    U = TIMEMATRIX(U,T);
end

patch = patches.patches;
interface = interfaces.interfaces;

n = numel(patches);
w = cell(1,n);
lambda = cell(1,n);
for k=1:n
    if isempty(patch{k}.timeSolver)
        w{k} = zeros(getnbddlfree(patch{k}.S),1);
        lambda{k} = zeros(getnbddlfree(interface{k}.S),1);
    else
        nbtimedof = getnbtimedof(patch{k}.timeSolver);
        sz_w = [getnbddlfree(patch{k}.S),nbtimedof];
        sz_lambda = [getnbddlfree(interface{k}.S),nbtimedof];
        w{k} = zeros(sz_w);
        w{k} = TIMEMATRIX(w{k},T);
        lambda{k} = zeros(sz_lambda);
        lambda{k} = TIMEMATRIX(lambda{k},T);
    end
end
% w = cellfun(@(patch) zeros(getnbddlfree(patch.S),1),patch,'UniformOutput',false);
% lambda = cellfun(@(interface) zeros(getnbddlfree(interface.S),1),interface,'UniformOutput',false);

% n = numel(patches);
% w = cell(1,n)
% lambda = cell(1,n);
% for k=1:n
%     w{k} = zeros(getnbddlfree(patch{k}.S),1);
%     lambda{k} = zeros(getnbddlfree(interface{k}.S),1);
% end

end
