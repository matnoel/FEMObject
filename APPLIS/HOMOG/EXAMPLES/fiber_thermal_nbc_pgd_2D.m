size_img=[1 1];
nb_elem=[20 30];

neg_side=10;
pos_side=1;
delta=0.01;
epsi=1;

[X,Y]=meshgrid(linspace(0,size_img(1),nb_elem(1)+1),linspace(0,size_img(2),nb_elem(2)+1));
r_fibre=0.3338;
ls=sqrt((X-0.5).^2+(Y-0.5).^2)-r_fibre;

homog=THERMALHOMOGNBCPGD(size_img, nb_elem, ls, neg_side, pos_side,delta,epsi,'tol',1e-5);
homog=homog.solve('residual',0,'maxorder',100,'update',0,'righthandSD',1);
homog=homog.findeffectiveconductivity('maxorder',100);

disp(homog.C)
homog=homog.computeflux('tol',1e-5);
homog.export('fiber_thermal_nbc_pgd_2D.vti')
