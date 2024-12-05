function av=variance(apc)

av=sum(double(apc.MYDOUBLE).^2,apc.pcdim);
av=av-mean(apc).^2;
