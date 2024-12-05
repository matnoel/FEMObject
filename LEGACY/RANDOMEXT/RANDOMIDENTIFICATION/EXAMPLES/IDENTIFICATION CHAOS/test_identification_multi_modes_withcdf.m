%% definition des lois bimodales
nbmodes=2;
x1 = setnumber(RVUNIFORM(0,nbmodes),1);
switch nbmodes
case 1

case 2 
    mu = [1,2];
    si = [0.2,0.4];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2);
    figure(10)
    clf


case 3
    mu = [0,1,2];
    si = [0.3,0.4,0.2];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2) + (x>=2 & x<=3).*mu(3);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2) + (x>=2 & x<=3).*si(3);
case 4
    mu = [-1,1,2,5];
    si = [0.2,0.4,0.1,0.2];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2) + (x>=2 & x<=3).*mu(3) + (x>=3 & x<=4).*mu(4);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2) + (x>=2 & x<=3).*si(3) + (x>=3 & x<=4).*si(4);
end

x2 = setnumber(CONDRANDVAR(@RVNORMAL,x1,@(x) mufun(x),@(x) sifun(x)),2);        

figure(10)
hold on
pdfplot(x2,'npts',100,'nbs',1e5,'ksdensity');

Xs = lhsrandom(x2,300,1)';
hold on
pdfsampleplot(Xs,'b','npts',70,'bar')



%%



