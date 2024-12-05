function f=calc_multivector(S,fun,varargin)
% function f=calc_vector(S,fun,varargin)
% S : MODEL
% fun : pointeur sur une methode de la classe ELEMENT
%       appel de  fe=fun(ELEMENT,S,varargin{:})
%       fe est un tableau de celulles
%       fe{i} est un MYDOUBLEND de taille nbddl-by-.-by-nbelem
%       ou nbddl est le nombre de ddl du groupe d'element i et nbelem son
%       nombre d'elements

display_ = ischarin('display',varargin);

if display_
    fprintf('\nCOMPUTING FINITE ELEMENT VECTOR\n')
end

if isa(fun,'char')
    fun = eval(['@' fun]);
end
fun = fcnchk(fun);

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);

for p=liste
    if display_
        fprintf('-> Computing element vectors of element group %3d / %3d ... ',p,S.nbgroupelem)
    end
    fe{p} = fun(S.groupelem{p},getnode(S),varargin{:});
    if display_
        fprintf('\n')
    end
end

if display_
    fprintf('-> Assembling multivector ... ')
end

f = assemble_multivectorelem(S,fe,varargin{:});

if display_
    fprintf('\n')
end
