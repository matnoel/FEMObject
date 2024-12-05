size_img=[1 1 1];
nb_elem=[50 50 5];
ls=LSCYLINDER(3,0.5,0.5,0,0,0,1,0.3338);
neg_side=1;
pos_side=1;

homog=THERMALHOMOGPBC(size_img, nb_elem, ls, neg_side, pos_side);
homog=homog.solve();
homog=homog.findeffectiveconductivity();
disp(homog.C)

homog.plot_ls()
homog.plot_mode_micro(1)
homog.plot_mode_macro(1)
homog.plot_mode_tot(1)
homog.plot_mode_micro(2)
homog.plot_mode_macro(2)
homog.plot_mode_tot(2)
homog.plot_mode_micro(3)
homog.plot_mode_macro(3)
homog.plot_mode_tot(3)
