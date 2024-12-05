function [u,v,a] = solveSystemCell(pb,varargin)
% function u = solveSystemCell(pb,varargin)
% function [u,v] = solveSystemCell(pb,varargin)
% function [u,v,a] = solveSystemCell(pb,varargin)
% Solves problem pb and returns solution u as a cell array if nargout==1
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
        u_val = getvalue(u);
        v_val = getvalue(v);
        if nargout==1
            u = cell(1,2);
            u{1} = u_val;
            u{2} = v_val;
        else
            u = u_val;
            v = v_val;
        end
    elseif pb.timeOrder==2
        if ~isfield(pb,'u0')
            pb.u0 = [];
        end
        if ~isfield(pb,'v0')
            pb.v0 = [];
        end
        [u,result,v,a] = ddsolve(pb.timeSolver,pb.b,pb.M,pb.A,[],pb.u0,pb.v0);
        u_val = getvalue(u);
        v_val = getvalue(v);
        a_val = getvalue(a);
        if nargout==1
            u = cell(1,3);
            u{1} = u_val;
            u{2} = v_val;
            u{3} = a_val;
        else
            u = u_val;
            v = v_val;
            a = a_val;
        end
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
    if nargout==1
        u = {u};
    end
end
