size_img=[1 1 1];
nb_elem=[10 20 30];

neg_side=10;
pos_side=1;
delta=0.01;
epsi=1;

[X,Y,Z]=meshgrid(linspace(0,size_img(1),nb_elem(1)+1), ... 
    linspace(0,size_img(2),nb_elem(2)+1), ...
    linspace(0,size_img(3),nb_elem(3)+1));
r_fibre=0.3338;
ls=sqrt((X-0.5).^2+(Y-0.5).^2)-r_fibre;

homog=THERMALHOMOGPBCPGD(size_img, nb_elem, ls, neg_side, pos_side, delta, epsi,'tol',1e-3);
homog=homog.solve('residual',0,'maxorder',50,'update',0);
homog=homog.findeffectiveconductivity('maxorder',50);
disp(homog.C)
homog=homog.computeflux('tol',1e-4);
homog=homog.computemicroflux('tol',1e-4);
homog.export('fiber_thermal_pbc_pgd_3D.vti');
