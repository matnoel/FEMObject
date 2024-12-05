function pb = calcOperator(pb,varargin)
% function pb = calcOperator(pb,varargin)
% Compute operator and update right-hand side of problem pb
% pb: struct containing mesh, operator, right-hand side, solver of the problem

if isfield(pb,'timeSolver') && ~isempty(pb.timeSolver)
    pb.M = calc_mass(pb.S);
end
if ~isfield(pb,'solver') || isempty(pb.solver)
    [pb.A,pb.b0] = calc_rigi(pb.S);
    pb.b0 = -pb.b0;
    if isfield(pb,'timeSolver') && ~isempty(pb.timeSolver)
        if isfield(pb,'loadFunction') && ~isempty(pb.loadFunction)
            pb.b0 = pb.b0*pb.loadFunction(pb.timeSolver);
        else
            pb.b0 = pb.b0*one(pb.timeSolver);
        end
    end
    if isfield(pb,'b') && ~isempty(pb.b)
        pb.b = pb.b + pb.b0;
    else
        pb.b = pb.b0;
    end
else
    if isa(pb.solver,'NEWTONSOLVER')
        pb.A = @(u) calc_fint(pb.S,u);
        pb.Atang = @(u) calc_rigitang(pb.S,u);
    end
end

end
