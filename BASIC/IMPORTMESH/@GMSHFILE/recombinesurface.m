function u = recombinesurface(u,number)
% function u = recombinesurface(u,number)

u = stepcounter(u);
s = ['Recombine Surface{' num2str(number) '} ;\n'];
u = addstring(u,s);
