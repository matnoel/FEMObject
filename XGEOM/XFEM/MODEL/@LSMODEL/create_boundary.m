function [frontiere,interfaceglob,interface] = create_boundary(M)
% function [frontiere,interfaceglob,interface] = create_boundary(M)

[frontiere,interfaceglob,interface] = create_boundary(M.MODEL);
frontiere = LSMODEL(frontiere,restrict(getlevelsets(M),M,frontiere));


