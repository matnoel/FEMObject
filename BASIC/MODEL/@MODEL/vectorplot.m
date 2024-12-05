function varargout = vectorplot(S,type,F,varargin)
% function varargout = vectorplot(S,type,u,,ampl)
% type : 'F' pour force , 'U' pour deplacement, ...
% u : vecteur des ddls 
% ampl : amplification

F = unfreevector(S,F);

ddl = DDL(DDLVECT(type,S.syscoord,'TRANS'));
switch type
    case 'F'
        name = 'Neumann';
        rep = findddl(S,ddl,'dual');
    case 'U'
        name = 'Dirichlet';
        rep = findddl(S,ddl);
    otherwise
        error('Wrong type')
end
if size(F,2)>1
    error('Wrong size of vector')
end

% if nargin>=4 && ~isempty(ampl)
%    F = ampl*F;
% end

F = reshape(F(rep,:),[length(DDL(ddl)),S.nbnode])';

H = quiver(FENODEFIELD(F),S,varargin{:});

% X = getcoord(S.node);
% X = mat2cell(X,size(X,1),ones(1,size(X,2)));
% F = mat2cell(F,size(F,1),ones(1,size(F,2)));
% 
% switch size(F,2)
%     case 2
%         H = quiver(X{:},F{:},varargin{:});
%     case 3
%         H = quiver3(X{:},F{:},varargin{:});
% end

ddlname = cell(1,length(ddl));
for k=1:length(ddl)
    if any(F(:,k))
        ddlname{k} = [getddlname(ddl,k), ','];
    end
end
ddlname = [ddlname{:}];
leg = [ddlname, name];

if nargout>=1
    varargout{1} = H;
end
if nargout>=2
    varargout{2} = leg;
end
