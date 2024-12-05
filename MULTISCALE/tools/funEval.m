function pb = funEval(pb,xi,varargin)
% function pb = funEval(pb,xi,varargin)
% Evaluation of problem pb at points xi
% pb: struct containing mesh, operator, right-hand side, solver of the problem
% xi: N-by-d double containing of the evaluations of random variables
% where
% N is the number of samples
% d is the parametric dimension (number of random variables)

if israndom(pb.S)
    pb.S = randomeval(pb.S,xi,varargin{:});
end

if isfield(pb,'A') && israndom(pb.A)
    pb.A = randomeval(pb.A,xi,varargin{:});
end

if isfield(pb,'b') && israndom(pb.b)
    pb.b = randomeval(pb.b,xi,varargin{:});
end

if isfield(pb,'Atang') && israndom(pb.Atang)
    pb.Atang = randomeval(pb.Atang,xi,varargin{:});
end

if isfield(pb,'M') && israndom(pb.M)
    pb.M = randomeval(pb.M,xi,varargin{:});
end

if isfield(pb,'C') && israndom(pb.C)
    pb.C = randomeval(pb.C,xi,varargin{:});
end

if isfield(pb,'b0') && israndom(pb.b0)
    pb.b0 = randomeval(pb.b0,xi,varargin{:});
end

if isfield(pb,'u0') && israndom(pb.u0)
    pb.b0 = randomeval(pb.u0,xi,varargin{:});
end

if isfield(pb,'v0') && israndom(pb.v0)
    pb.b0 = randomeval(pb.v0,xi,varargin{:});
end

end
