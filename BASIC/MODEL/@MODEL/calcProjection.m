function [P,numnode2] = calcProjection(M1,M2,ddls,varargin)
% function [P,numnode2] = calcProjection(M1,M2,ddls,varargin)
% Calculates projection operator from model M1 to model M2
% M1, M2 : MODEL (M1: fine mesh, M2: coarse mesh)
% ddls : liste des ddls a bloquer {'UX','U','R',...} pour la mecanique,
%        'T' pour la thermique ou autre ... selon le choix du MODEL 
% ddls = 'all' par defaut
% P : Projection operator

if nargin<3 || isempty(ddls) || strcmp(ddls,'all')
    ddls = 'all';
elseif isa(ddls,'char')
    ddls = {ddls};
end
full = getcharin('full',varargin,false);
free = getcharin('free',varargin,true);

I = speye(M2.nbddl,M2.nbddl);

if full
    numnode2 = 1:M2.nbnode;
else
    [~,numnode2] = intersect(M2,M1,'strict',0);
end
numddl2 = findddl(M2,ddls,numnode2);
ddl = getddl(M2.node,getnumber(M2.node));
Ptmp = eval_sol(M2,I(:,numddl2),M1,ddl);
% P = sparse(squeeze(Ptmp));
Ptmp  = permute(Ptmp,[2,1,3]);
Ptmp = Ptmp(:,:);
P = sparse(M2.nbddl,M1.nbddl);
P(numddl2,:) = Ptmp;

if free
    P = freevector(M1,P,2);
    P = freevector(M2,P,1);
end

end
