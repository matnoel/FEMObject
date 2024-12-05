function A = calcIterationOperator(W,glob,patches,interfaces,V,w)
% function A = calcIterationOperator(W,glob,patches,interfaces,V,w)
% Calculates iteration operator A(W) = W - Phi(Psi(W),W)
% W: a global solution
% glob: Glob
% patches: Patches
% interfaces: Interfaces
% V: current global solution V_{k}
% w: previous local solution w_{k-1}
% vV: current global velocity vV_{k}
% aV: current global acceleration aV_{k}

if nargin<5 || (isempty(V) && isempty(w))
    [V,w,~] = initializeVariables(glob,patches,interfaces);
elseif isempty(V)
    [V,~,~] = initializeVariables(glob,patches,interfaces);
elseif nargin<6 || isempty(w)
    [~,w,~] = initializeVariables(glob,patches,interfaces);
end

if ~isempty(glob.timeSolver)
    T = gettimemodel(glob.timeSolver);
    sz = [getnbddlfree(glob.S),getnbtimedof(glob.timeSolver)];
    W = reshape(W,sz);
    W = TIMEMATRIX(W,T);
end

% Local problems
% Local solutions (Theta{k},Psi{k}) without change of variable
%                           Zeta{k} with change of variable
n = numel(patches);
patch = patches.patches;
interface = interfaces.interfaces;
Psi = cell(1,n);
for k=1:n
    patch{k}.b = zeros(size(patch{k}.b));
    if ~isempty(patch{k}.timeSolver)
        T = gettimemodel(patch{k}.timeSolver);
        patch{k}.b = patch{k}.b*zero(T);
    end
    patch{k}.initializationType = 'zero';
    patch{k} = patch{k}.linearizeOperator(w{k});
    if isempty(patch{k}.timeSolver)
        [~,Psi{k}] = patch{k}.solve(interface{k},W);
    else
        if patch{k}.timeOrder==1
            [~,Psi{k}] = patch{k}.solve(interface{k},W);%,[],[],vW,[],[]);
        elseif patch{k}.timeOrder==2
            [~,Psi{k}] = patch{k}.solve(interface{k},W);%,[],[],vW,[],[],aW,[],[]);
        end
    end
end

% Global problem
% Global solution Phi
g = glob;
g.b_out = zeros(size(g.b_out));
if ~isempty(g.timeSolver)
    T = gettimemodel(g.timeSolver);
    g.b_out = g.b_out*zero(T);
end
g.initializationType = 'zero';
g = g.linearizeOperator(V);
if isempty(g.timeSolver)
    Phi = g.solve(interfaces,Psi,W);
else
    if g.timeOrder>=1
        vW = diff(g.timeSolver,W);
    end
    if g.timeOrder>=2
        aW = diff(g.timeSolver,vW);
    end
    if isadsolver(g.timeSolver) && g.timeOrder==1
        Phi = glob.solve(interfaces,Psi,W,[],vW,[]);
    elseif isaddsolver(g.timeSolver) && g.timeOrder==2
        Phi = glob.solve(interfaces,Psi,W,[],vW,[],aW,[]);
    end
end

% Iteration operator A
A = W - Phi;
if ~isempty(glob.timeSolver)
    A = getvalue(A);
    A = reshape(A,prod(sz),1);
end

end
