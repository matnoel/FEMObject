function plottext(M,varargin)
xplot=double(getcoord(M.node));
xplot=POINT(mean(xplot,1));
plottext(xplot,varargin{:}); 

