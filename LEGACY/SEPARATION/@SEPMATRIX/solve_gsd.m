function [u,varargout] = solve_gsd(dimgsd,A,b,solver,solverlambda,solverlambdaupdate)
% function [u,varargout] = solve_gsd(dimgsd,A,b,solver,solverlambda,solverlambdaupdate)

varargout=cell(1,nargout-1);

if nargin<6
    solverlambdaupdate=solverlambda;
end

if dimgsd~=1
    dim = (1:getdim(A));
    dim(1)=dimgsd;
    dim(dimgsd)=1;
    A = permutedim(A,dim);
    b = permutedim(b,dim);
    if strcmp(getparam(solver,'errorindicator'),'reference')
        solref = getparam(solver,'reference');
        solref = permute(solref,dim);
        solver = setparam(solver,'reference',solref);
    end
    
    [u,varargout{:}] = solve_gsd(1,A,b,solver,solverlambda,solverlambdaupdate);
    u = permutedim(u,dim);
    return
end

switch getparam(solver,'type')
    case 'alterne'
        [u,varargout{:}] = solve_gsd_alterne(A,b,solver,solverlambda,solverlambdaupdate);
    case 'arnoldi'
        [u,varargout{:}] = solve_gsd_arnoldi(A,b,solver,solverlambda,solverlambdaupdate);
    otherwise
        error('mauvais type de solver')
end
