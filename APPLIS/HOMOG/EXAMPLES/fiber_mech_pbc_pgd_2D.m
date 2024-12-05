size_img=[1 1];
nb_elem=[50 60];
delta=0.01;

[X,Y]=meshgrid(linspace(0,size_img(1),nb_elem(1)+1),linspace(0,size_img(2),nb_elem(2)+1));
r_fibre=0.3338;
lset=sqrt((X-0.5).^2+(Y-0.5).^2)-r_fibre;

%%
E_ns=10;
E_ps=1;
nu_ns=0.3;
nu_ps=0.3;

C_ns=elasticity_tensor(E_ns,nu_ns,2);
C_ps=elasticity_tensor(E_ps,nu_ps,2);
clear homog;
epsilon=1;
homog=MECHHOMOGPBCPGD(size_img, nb_elem, lset, C_ns, C_ps,delta,epsilon,'tol',1e-6,'maxorder',100);
clear lset;
%%
homog=homog.solve('residual',0,'maxorder',200,'update',0,'tol',2e-4);

%%
homog=homog.computefields('update',1,'updatestep',25,'tol',1e-6);

%%
homog=homog.findeffectiveproperties();
disp(homog.C)

homog.export('fiber_mech_pbc_pgd_2D.vti');
