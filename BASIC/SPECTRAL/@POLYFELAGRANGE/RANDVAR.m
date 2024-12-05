function rv=RANDVAR(h)
domain = get(h,'domain');
rv=RVUNIFORM(domain(1),domain(2));
rv = setnumber(rv,getnumber(h));