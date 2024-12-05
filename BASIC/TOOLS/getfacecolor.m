function c = getfacecolor(i)

k=floor(i/7);
j=mod(i-1,7)+1;
C = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;1 1 1];
c = C(j,:)/(k+1);

