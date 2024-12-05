function G = gmshfile(L,cl,numberpoints,numberline)
% function G = gmshfile(L,cl,numberpoints,numberline)
% L : LIGNE
% cl : characteristic length

if nargin<=2
    numberpoints = 1:2;
    numberline = 1;
end

G = GMSHFILE();
P = getvertices(L);
G = createpoints(G,P(1:2),cl,numberpoints);
G = createline(G,numberpoints,numberline);
