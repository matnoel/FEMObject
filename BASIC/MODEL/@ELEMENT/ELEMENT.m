function elem = ELEMENT(varargin)
% function elem = ELEMENT()
%
% ddlnode  : nom des ddl aux noeuds (exemple : {'UX','UY'} )
% ddlnodedual  : nom des ddl duaux aux noeuds (exemple : {'FX','FY'} )
% ddlgauss : nom des ddl aux points de gauss (exemple : {'SMXX','SMYY','SMXY'} )
% ddlgaussprimal : nom des ddl aux points de gauss (exemple : {'EPXX','EPYY','EPXY'} )
%
% function elem = ELEMENT('material',mat,'option',option)
% mat : MATERIAL
% option : 'DEFO' , 'CONT' , 'BORD'

mat = getclassin('MATERIAL',varargin);
%if isa(mat,'MATERIAL') & isempty(getnumber(mat))
%    error('le materiau doit avoir un numero')
%end
option = getcharin('option',varargin,' ');
parent = getcharin('parent',varargin,' ');

elem.param = struct();
elem.parent = parent;
elem.local = 0;

elem.material=mat;
elem.option=option;
elem.ddlnode=DDL();
elem.nbddlpernode=0;
elem.ddlnodedual=DDL();
elem.ddlgauss=DDL();
elem.nbddlpergauss=0;
elem.ddlgaussdual=DDL();
elem.nbddl=0;
elem.numddl=[];

elem.lstype = 'indomain';
elem.lsenrich = 0 ;
elem.lsnumber = [];
elem.lsnature = [];

elem.param = setfield(elem.param,'initializeBN',false);

elem=class(elem,'ELEMENT');



