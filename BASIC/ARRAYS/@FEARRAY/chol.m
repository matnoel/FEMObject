function R = chol(K)

rep = K.ddlfree;
repb = K.ddlbloque;

R = zeros(size(K));
R(rep,rep) = chol(double(K.MYDOUBLE(rep,rep)));

R = FEMATRIX(R,repb);
