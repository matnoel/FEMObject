function [U,w,lambda] = initializeRandomVariables(glob,patches,interfaces,bases)
% function [U,w,lambda] = initializeRandomVariables(glob,patches,interfaces,bases)
% Initializes global solution U, local solution w and Lagrange multipier lambda
%
% Inputs:
% glob: Global
% patches: Patches
% interfaces: Interfaces
% bases: FunctionalBases
%
% Outputs:
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

d = length(bases);
I = MultiIndices(zeros(1,d));
basis = SparseTensorProductFunctionalBasis(bases,I);

if isempty(glob.timeSolver)
    sz = getnbddlfree(glob.S);
else
    sz = [getnbddlfree(glob.S),getnbtimedof(glob.timeSolver)];
end
data = zeros([cardinal(basis),sz]);
U = FunctionalBasisArray(data,basis,sz);

patch = patches.patches;
interface = interfaces.interfaces;

n = numel(patches);
sz = cell(1,n);
for k=1:n
    if isempty(patch{k}.timeSolver)
        sz{k} = getnbddl(patch{k}.S);
    else
        sz{k} = [getnbddl(patch{k}.S),getnbtimedof(patch{k}.timeSolver)];
    end
end
% sz = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
data = cellfun(@(sz) zeros([cardinal(basis),sz]),sz,'UniformOutput',false);
w = cellfun(@(data,sz) FunctionalBasisArray(data,basis,sz),data,sz,'UniformOutput',false);

for k=1:n
    if isempty(patch{k}.timeSolver)
        sz{k} = getnbddl(interface{k}.S);
    else
        sz{k} = [getnbddl(interface{k}.S),getnbtimedof(patch{k}.timeSolver)];
    end
end
% sz = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
data = cellfun(@(sz) zeros([cardinal(basis),sz]),sz,'UniformOutput',false);
lambda = cellfun(@(data,sz) FunctionalBasisArray(data,basis,sz),data,sz,'UniformOutput',false);

% n = numel(patches);
% w = cell(1,n);
% lambda = cell(1,n);
% for k=1:n
%     if isempty(patch{k}.timeSolver)
%         sz = getnbddl(patch{k}.S);
%     else
%         sz = [getnbddl(patch{k}.S),getnbtimedof(patch{k}.timeSolver)];
%     end
%     data = zeros(cardinal(basis),sz);
%     w{k} = FunctionalBasisArray(data,basis,sz);
%     if isempty(patch{k}.timeSolver)
%         sz = getnbddl(interface{k}.S);
%     else
%         sz = [getnbddl(interface{k}.S),getnbtimedof(patch{k}.timeSolver)];
%     end
%     data = zeros(cardinal(basis),sz);
%     lambda{k} = FunctionalBasisArray(data,basis,sz);
% end

end
