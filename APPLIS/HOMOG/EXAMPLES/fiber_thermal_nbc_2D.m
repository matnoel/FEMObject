%See Jiang/Ostoja 2002

size_img=[1 1];
nb_elem=[100 100];
ls=LSCIRCLE(0.5,0.5,0.3338);
neg_side=10;
pos_side=1;

homog=THERMALHOMOGNBC(size_img, nb_elem, ls, neg_side, pos_side);
homog=homog.solve();
homog=homog.findeffectiveconductivity();
disp(homog.C)

homog.plot_ls()
homog.plot_mode(1)
homog.plot_mode(2)
