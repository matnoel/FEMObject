function u0=calc_init_dirichlet(BC)

u0 = zeros(getnbddlfree(BC),1);
u0 = unfreevector(BC,u0);