function M = calc_matrix(S,fun,varargin)
% function M = calc_matrix(S,fun,varargin)
% S : MODEL
% fun : pointeur sur une methode de la classe ELEMENT
%       appel de  me=fun(ELEMENT,S,varargin{:})
%       me est un tableau de celulles
%       me{i} est un MYDOUBLEND de taille nbddl-by-nbddl-by-nbelem
%       ou nbddl est le nombre de ddl du groupe d'element i et nbelem son
%       nombre d'elements

display_ = ischarin('display',varargin);
if display_
    fprintf('\n COMPUTING FINITE ELEMENT MATRIX \n')
end

if isa(fun,'char')
    fun = eval(['@' fun]);
end
fun = fcnchk(fun);

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);

me = cell(1,getnbgroupelem(S));
tic
for p=liste
    if display_
        fprintf('-> Computing element matrices of element group %3d / %3d ... ',p,S.nbgroupelem)
    end
    me{p} = fun(S.groupelem{p},getnode(S),varargin{:});
    if display_
        fprintf('Elapsed time is %.3f seconds.',toc)
        fprintf('\n')
    end
    
end

if display_
    fprintf('-> Assembling matrix ... ')
end
M = assemble_matrixelem(S,me,varargin{:});
if display_
    fprintf('Elapsed time is %.3f seconds.\n',toc)
    fprintf('  Total elapsed time is %.3f seconds.\n',toc)
end

