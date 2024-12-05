function [M]=calc_pcmatrix(S,PC,fun,varargin)
% function M=calc_multimatrix(S,fun,varargin)
% S : MODEL
% fun : pointeur sur une methode de la classe ELEMENT
%       appel de  me=fun(ELEMENT,S,varargin{:})
%       me est un tableau de celulles
%       me{i} est la matrice elementaire de taille ELEMENT.nbddl * ELEMENT.nbddl
%       associee a l'element i de ELEMENT
% nmat : nombre de matrices calculees simultanement

display_ = ischarin('display',varargin);

if display_
    fprintf('\nCOMPUTING FINITE ELEMENT MATRIX\n')
end

if isa(fun,'char')
    fun = eval(['@' fun]);
end
fun = fcnchk(fun);

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);

tic
for p=liste
    if display_
        fprintf('\n-> Computing element matrices of element group %3d / %3d ... ',p,S.nbgroupelem)
    end
    me{p} = fun(S.groupelem{p},getnode(S),PC,varargin{:});
    if display_
        fprintf('Elapsed time is %.3f seconds.',toc)
    end
end

M = assemble_pcmatrixelem(S,me,'selgroup',liste);

if ~israndom(M)
    M = PCRADIALMATRIX({M},size(M),one(PC));
end

if display_
    fprintf('Total elapsed time is %.3f seconds.\n',toc)
end
