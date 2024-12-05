function A = calcIterationOperatorRandom(W,glob,patches,interfaces,s,bases,ls,rv,V,w)
% function A = calcIterationOperatorRandom(W,glob,patches,interfaces,s,bases,ls,rv,V,w)
% Calculates iteration operator A(W) = W - Phi(Psi(W),W)
% W: FunctionalBasisArray of global solution W
% glob: Glob
% patches: Patches
% interfaces: Interfaces
% s: AdaptiveSparseTensorAlgorithm
% bases: FunctionalBases (should provide orthonormal bases) or
% FullTensorProductFunctionalBasis
% ls: LinearModelLearningSquareLoss
% rv: RandomVector or RandomVariable (optional)
% V: FunctionalBasisArray of current global solution V_{k}
% w: FunctionalBasisArray of current local solution w_{k-1}

if isa(bases,'FullTensorProductFunctionalBasis')
    bases = bases.bases;
elseif isa(bases,'FunctionalBasis')
    bases = FunctionalBases({bases});
end

d = length(bases);
rvb = getRandomVector(bases);
if nargin<7 || isempty(rv)
    rv = rvb;
elseif isa(rv,'RandomVariable')
    rv = RandomVector(rv,d);
end
if nargin<8 || (isempty(V) && isempty(w))
    [V,w,~] = initializeRandomVariables(glob,patches,interfaces,bases);
elseif isempty(V)
    [V,~,~] = initializeRandomVariables(glob,patches,interfaces,bases);
elseif nargin<9 || isempty(w)
    [~,w,~] = initializeRandomVariables(glob,patches,interfaces,bases);
end

n = numel(patches);

indices = V.basis.indices;
for k=1:n
    indices = indices.addIndices(w{k}.basis.indices);
end
basis = SparseTensorProductFunctionalBasis(bases,indices);
sz = getnbddlfree(glob.S);
W = FunctionalBasisArray(W,basis,sz);

% Local problems
% Local solutions (Theta{k},Psi{k}) without change of variable
%                           Zeta{k} with change of variable
Psi = cell(1,n);
% err_Psi = cell(1,n);
% dim_Psi = cell(1,n);
% N = cell(1,n);
% y_Psi = cell(1,n);
s.display = false;
s.displayIterations = false;
ls.errorEstimation = true;
patch = patches.patches;
interface = interfaces.interfaces;
parfor k=1:n
    patchk = patch{k};
    patchk.b = zeros(size(patchk.b));
    patchk.initializationType = 'zero';
    w_proj = w{k}.projection(basis);
    
    fun = @(xi) solve(linearizeOperator(calcOperator(patchk.eval(xi)),...
        w_proj(xi)'),interface{k},W(xi)');
    fun = UserDefinedFunction(fun,d,2);
    fun.evaluationAtMultiplePoints = false;
    
    if s.display || s.displayIterations
        fprintf('\n+-----------+------------+------------+\n');
        fprintf('| Dim basis | Nb samples |  CV error  |\n');
        fprintf('+-----------+------------+------------+\n');
    end
    
    m = fun.outputSize;
    H = repmat({basis},m,1);
    % H = cell(m,1); [H{:}] = deal(basis);
    % H = cell(m,1); H(:) = {basis};
    N = cardinal(basis);
    x = random(rv,N);
    fun = @(x) fun.evalCell(x);
    y = fun(x);
    if eq(rv,rvb)
        xb = x;
    else
        xb = transfer(rv,rvb,x);
    end
    A = cellfun(@(basis) basis.eval(xb),H,'UniformOutput',false);
    
    u = cell(m,1);
    err = repmat({Inf},m,1);
    [u,err,x,y] = s.adaptSamplingCell(fun,H,ls,rv,u,err,x,y,A);
    
    if s.display && ~s.displayIterations
        for i=1:m
            fprintf('| %9d | %10d | %4.4e |\n',cardinal(H{i}),size(x,1),norm(err{i}));
        end
    end
    
    if s.display || s.displayIterations
        fprintf('+-----------+------------+------------+\n');
    end
    
    if s.fullOutput
        u = fullOutputConversion(s,u);
    end
    
    % [u,err,x,y] = s.leastSquaresCell(fun,bases,ls,rv);
    
    [~,Psi{k}] = u{:};
    % [~,err_Psi{k}] = err{:};
    % [~,y_Psi{k}] = y{:};
    % dim_Psi{k} = cardinal(Psi{k}.basis);
    % N{k} = size(x,1);
end

% Global problem
% Global solution Phi
indices = W.basis.indices;
for k =1:n
    indices = indices.addIndices(Psi{k}.basis.indices);
end
basis = SparseTensorProductFunctionalBasis(bases,indices);
W_proj = W.projection(basis);
Psi_proj = cellfun(@(x) x.projection(basis),Psi,'UniformOutput',false);
g = glob;
g.b_out = zeros(size(g.b_out));
g.initializationType = 'zero';
g = g.linearizeOperator(V);
Phi = g.solve(interfaces,cellfun(@(x) x.data',Psi_proj,'UniformOutput',false),...
    W_proj.data')';
Phi = FunctionalBasisArray(Phi,basis,sz);

% Iteration operator A
A = W_proj.data - Phi.data;
A = double(A);
A = A(:);

end
