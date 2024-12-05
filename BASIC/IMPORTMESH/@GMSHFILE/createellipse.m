function u = createellipse(u,numcenter,numpoints,numbermajorpoint,numberline)
% function u = createellipse(u,numcenter,numpoints,numbermajorpoint,numberline)

if length(numpoints)~=2
    error('pour creer un arc d''elipse, il faut 2 numeros de points')
end
u = createentity(u,'Ellipse',[numpoints(1), numcenter, numbermajorpoint , numpoints(2)],numberline);
