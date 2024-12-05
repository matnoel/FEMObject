function [u,v,a] = solveSystem(pb,varargin)
% function u = solveSystem(pb,varargin)
% function [u,v] = solveSystem(pb,varargin)
% function [u,v,a] = solveSystem(pb,varargin)
% Solves problem pb
% pb: struct containing mesh, operator, right-hand side, solver of the problem
% u: solution
% For first- and second-order time-dependent problems,
% v: velocity
% For second-order time-dependent problems,
% a: acceleration

if isfield(pb,'timeSolver') && ~isempty(pb.timeSolver) && isfield(pb,'timeOrder') && ~isempty(pb.timeOrder)
    if pb.timeOrder==1
        if ~isfield(pb,'u0')
            pb.u0 = [];
        end
        [u,result,v] = dsolve(pb.timeSolver,pb.b,pb.M,pb.A,pb.u0);
        u = getvalue(u);
        v = getvalue(v);
    elseif pb.timeOrder==2
        if ~isfield(pb,'u0')
            pb.u0 = [];
        end
        if ~isfield(pb,'v0')
            pb.v0 = [];
        end
        [u,result,v,a] = ddsolve(pb.timeSolver,pb.b,pb.M,pb.A,[],pb.u0,pb.v0);
        u = getvalue(u);
        v = getvalue(v);
        a = getvalue(a);
    end
else
    if ~isfield(pb,'solver') || isempty(pb.solver)
        u = pb.A\pb.b;
    else
        if isa(pb.solver,'NEWTONSOLVER')
            inittype = getcharin('inittype',varargin,'zero');
            u0 = setInitialSolution(inittype,pb.b);
            % fsolver = @(A,b) solve(A,b);
            % u = solve(pb.solver,pb.b,pb.A,pb.Atang,u0,[],fsolver);
            u = solve(pb.solver,pb.b,pb.A,pb.Atang,u0);
        end
    end
end
