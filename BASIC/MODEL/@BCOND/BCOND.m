function BC = BCOND(varargin)
% constructeur de la classe BCOND
% function BC = BCOND(node)

if ~isa(varargin{1},'NODE')
    error('entrer un objet NODE')
end

BC.nbddl = getnbddl(varargin{1});

BC.nbcl=0;
BC.BC=cell(0,1);
BC.ddlbloque = [];
BC.ddlfree = [1:BC.nbddl];
BC.ddlnoenrich = [1:BC.nbddl];
BC.ddlenrich = [];

groupddl = getgroupddl(varargin{1});
for p=1:length(groupddl)
    rep = findenrich(groupddl{p}.ddlnode);
    BC.ddlenrich = union(BC.ddlenrich,unique(groupddl{p}.numddl(:,rep)));
end
BC.ddlnoenrich = setdiff(BC.ddlnoenrich,BC.ddlenrich);

BC = class(BC,'BCOND');
