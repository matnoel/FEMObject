function [discon2con,con2discon] = calcTransferDiscontinuous(model)
% [discon2con,con2discon] = calcTransferDiscontinuous(model)
%_discon2con: sparse operator such that discon2con*discontinuous=continuous
%_con2discon: index such that continuous(con2discon)=discontinuous

tolerance = getfemobjectoptions('tolerancepoint') ;

% Get coordinates of continuous mesh
coordDiscontinuous = getDomainCoord(model) ;
coordDiscontinuous = round(coordDiscontinuous/tolerance)*tolerance ;
[coordContinuous,~,con2discon] = unique(coordDiscontinuous,'rows','stable') ;

% Compare: coordContinuous(matchedC,:)=coordDiscontinuous(d2c,:)
[matchedC,d2c] = ismember(coordContinuous,coordDiscontinuous,'rows') ;
if ~all(matchedC)
    warning('Not all continuous coordinates were matched. Tolerance too low?');
end

% Nodes not listed in d2c are duplicates. 
% List them and find matching nodes in continuous mesh.
duplicate = setdiff((1:size(coordDiscontinuous,1))',d2c) ;
[matchedD,c2d] = ismember(coordDiscontinuous,coordContinuous,'rows') ;
if any(~matchedD)
    warning('Not all discontinuous coordinates were matched. Tolerance too low?');
end

% Build matrix such that d2cOp(m,n) = 1 if
% coordDiscontinuous(n,:)==coordContinuous(m,:) and else 0.
i = [find(matchedC) ; c2d(duplicate)] ;
j = [d2c ; duplicate] ;
discon2con = sparse(i,j,1,size(coordContinuous,1),size(coordDiscontinuous,1)) ;

% Normalise so that the sum of each line be one (i.e. smooth values by averaging)
discon2con = diag(sum(discon2con,2).^-1)*discon2con ;
end