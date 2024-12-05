function varargout = plot(L,varargin)
% function varargout = plot(L,varargin)

P1 = getcoord(L.P{1});
P2 = getcoord(L.P{2});
nodecoord = double([P1;P2]) ;
n = size(nodecoord);
nodecoord = [nodecoord,zeros(n(1),n(2)==1)];
connec = 1:2;
edgecolor = getcharin('color',varargin,'k') ;

H = patch('faces',connec,'vertices',nodecoord,'edgecolor',edgecolor);

if nargout>=1
    varargout{1} = H;
end
  