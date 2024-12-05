function u = NODE(varargin)
% function u = NODE(POINT,number)
% constructeur d'une liste de noeuds
% utiliser addnode pour ajouter des noeuds a cette liste

if nargin==0
    u = NODE(POINT());
elseif nargin==1 && isa(varargin{1},'NODE')
    u=varargin{1};
elseif nargin==1 &&  isa(varargin{1},'SYSCOORD')
    u.number = zeros(0,1);
    u.nbnode = 0 ;
    u.groupddl=cell(0,1);
    u.nbgroupddl = 0;
    u.repnodeingroupddl = zeros(u.nbnode,2);
    u.lsenrichnode = [];
    u.lsenrichtype = [];
    u.lsnumber = [];
    u.lsenrichnature = cell(0,1);
    u.repnature = [];    
    u.nbddl =0 ;
    u=class(u,'NODE',POINT(zeros(1,getindim(varargin{1}),0)));
    superiorto('POINT');
    inferiorto('ELEMENT');
elseif nargin==1 && (isa(varargin{1},'POINT') || isa(varargin{1},'double'))
    if isa(varargin{1},'POINT')
        u = NODE(varargin{1},[1:numel(varargin{1})]');
    elseif isa(varargin{1},'double')
        u = NODE(varargin{1},[1:size(varargin{1},1)]');    
    end
elseif nargin==2 
    P = POINT(varargin{1});
    u.number = reshape(varargin{2},numel(varargin{2}),1);
    u.nbnode = numel(u.number);
    u.groupddl=cell(0,1);
    u.nbgroupddl = 0;
    u.repnodeingroupddl = zeros(u.nbnode,2);
    u.lsenrichnode = [];
    u.lsenrichtype = [];
    u.lsnumber = [];
    u.lsenrichnature = cell(0,1);
    u.repnature = [];
    u.nbddl =0 ;
    u=class(u,'NODE',P);
    superiorto('POINT');
    inferiorto('ELEMENT');
else
    error('mauvais arguments pour NODE')
end
