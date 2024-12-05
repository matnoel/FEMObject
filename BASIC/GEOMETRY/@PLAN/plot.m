function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P1 = double(getcoord(D.P{1}));
P2 = double(getcoord(D.P{2}));
P3 = double(getcoord(D.P{3}));
nodecoord = [P1;P2;P3] ;
connec = [1,2,3];
edgecolor = getcharin('color',varargin,'k') ;

H = patch('faces',connec,'vertices',nodecoord,'edgecolor',edgecolor);
  
if nargout>=1
    varargout{1} = H;
end
