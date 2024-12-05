
function varargout=contourplot(c,D,varargin)

if israndom(c)
    warning('la levelset est aleatoire : pas de contourplot possible')
    return
end
contourval=getcharin('value',varargin,0);
c = lseval(c,D);
D = setlevelsets(D,c);
D = lssplitelem(D);


for i=1:length(contourval)
    
    for p=1:getnbgroupelem(D)
        elem = getgroupelem(D,p);
        lstype = getlstype(elem);
        if strcmp(lstype,'cut') || strcmp(lstype,'bicut') || isenrich(elem)
            Dadd = MODEL(getmode(D));
            
            switch lstype
                case {'cut'}
                    [elem,node]=contour(elem,D.node,getvalue(getlssupport(c)),contourval(i));
                    Dadd = addelem(Dadd,elem);
                    Dadd = addnode(Dadd,node);
                    %plot(Dadd,'facecolor','g','facelighting','gouraud','facealpha',.1,'color','k')
                    %plot(Dadd,'facecolor','none','facelighting','gouraud','facealpha',.1,'color','k')
                    
                    % keyboard
                case 'bicut'
                    
                    lstip = getlstip(c,getparam(getgroupelem(D,p),'tipnumber'));
                    lssupport = getlssupport(c);
                    [elemcutin,elemcutout,nodeplus,xnodein,xnodeout,ls1value]=...
                        lsdivideelem(elem,lstip,D.node,getvalue(lssupport));
                    
                    %figure
                    %clf
                    %plot(D,'selgroup',3,'facecolor','none','facelighting','gouraud','facealpha',.5,'color','k')
                    %plot(D,'selgroup',1,'facecolor','r','facelighting','gouraud','facealpha',.5,'color','k')
                    %plot(D,'selgroup',2,'facecolor','g','facelighting','gouraud','facealpha',.5,'color','k')
                    %keyboard
                    
                    for k=1:length(elemcutin)
                        [elem,node]=contour(elemcutin{k},nodeplus,ls1value,contourval(i));
                        Dadd = addelem(Dadd,elem);
                        Dadd = addnode(Dadd,node);
                    end
                    
                    %plot(Dadd,'facecolor','none','facelighting','gouraud','facealpha',.5,'color','b')
                    % Dadd = MODEL(getmode(D));
                    %  for k=1:length(elemcutout)
                    % [elem,node]=contour(elemcutout{k},nodeplus,ls1value,contourval(i));
                    % Dadd = addelem(Dadd,elem);
                    % Dadd = addnode(Dadd,node);
                    % end
                    % plot(Dadd,'facecolor','w','facelighting','gouraud','facealpha',.5,'color','m')
                    
                    %keyboard
                    
                    
                    
                otherwise
                    %warning('mal programme : prolonge la fissure')
                    %[elem,node]=contour(elem,D.node,getvalue(c.LEVELSETS{1}),contourval(i));
                    %   Dadd = addelem(Dadd,elem);
                    %   Dadd = addnode(Dadd,node);
                    
            end
            
            if length(contourval)==1 | ischarin('color',varargin)
                col = getcharin('color',varargin,'r');
                varargin = setcharin('color',varargin,col);
            else
                col = contourval(i);
            end
            plot(Dadd,'color',col,varargin{:});
            
        end
    end
end

if nargout>=1
    varargout{1}=node;
    varargout{2}=coseg;
end

%
%
% S.ls = LEVELSETS(c);
% S = lssplitelem(S,'lsenrichtype',2);
% S1 = keepgroupelem(S,'cut');
% contourplot(getlevelset(c,1),S1,varargin{:});
%
% S2 = keepgroupelem(S,'bicut');
%
% if getnbelem(S2)>0
%     warning('fissure mal representee en pointe')
%     contourplot(getlevelset(c,1),S2,varargin{:});
% end
