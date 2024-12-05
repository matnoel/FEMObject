function elem = ELEMENTGEOM(varargin)
% function elem = ELEMENTGEOM(dim,node,numelem,connec,syscoordlocal,syscoord)
% dim : dimension de l'element
% node     : objet NODE ;
% numelem  : numero des element
% connec     : table de connectivite ;
% param : parametres des elements
if nargin==0
    elem=struct('dim',[],'numelem',[],'connec',[],'nbelem',[],...
        'nbnode',[],'syscoordlocal',SYSCOORD(),'syscoord',[]);
    elem=class(elem,'ELEMENTGEOM');

elseif nargin==1 && isa(varargin{1},'ELEMENT')
    elem=getelemgeom(varargin{1});
elseif nargin==1 && isa(varargin{1},'ELEMENTGEOM')
    elem=varargin{1};
elseif nargin==1 && isa(varargin{1},'double')
    elem=struct('dim',varargin{1},'numelem',[],'connec',[],'nbelem',[],...
        'nbnode',[],'syscoordlocal',SYSCOORD(),'syscoord',[]);
    elem=class(elem,'ELEMENTGEOM');    
else
    elem.dim=varargin{1};
    elem.numelem = varargin{3}(:);
    elem.connec=varargin{4};
    elem.nbelem = size(elem.connec,1);
    elem.nbnode=size(elem.connec,2);
    elem.syscoordlocal=varargin{5};
    if nargin<6 || isempty(varargin{6})
        elem.syscoord = getsyscoord(varargin{2});
    else
        elem.syscoord = varargin{6};
    end

    elem=class(elem,'ELEMENTGEOM');

end
