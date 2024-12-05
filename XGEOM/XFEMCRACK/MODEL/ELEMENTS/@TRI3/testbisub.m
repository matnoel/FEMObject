function testbisub(t,ls1,ls2)

[cinin,cinout,coutin,coutout,xlnodeplus]=bilsdivide_oneelem(ls1,ls2);%
nodecoord = [nodelocalcoord(t);xlnodeplus];

t = TRI3(NODE(nodelocalcoord(t)));
figure(1)
clf
subplot(3,2,1)
patch('faces',cinin,'vertices',nodecoord,'facecolor','r','edgecolor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('inin')
subplot(3,2,2)
patch('faces',cinout,'vertices',nodecoord,'facecolor','y','edgecolor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('inout')
subplot(3,2,3)
patch('faces',coutin,'vertices',nodecoord,'facecolor','g','edgecolor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('outin')

subplot(3,2,4)
patch('faces',coutout,'vertices',nodecoord,'facecolor','m','edgecolor','k');
plot(t,NODE(nodelocalcoord(t)))
axis square
title('outout')

[ginin,ginout,goutin,goutout]=calc_bilssubgauss(t,ls1,ls2,1);

subplot(3,2,5:6)
patch('faces',cinin,'vertices',nodecoord,'facecolor','r','edgecolor','k');
plot(POINT(ginin.coord),'yo');
patch('faces',cinout,'vertices',nodecoord,'facecolor','y','edgecolor','k');
plot(POINT(ginout.coord),'ro');
patch('faces',coutin,'vertices',nodecoord,'facecolor','g','edgecolor','k');
plot(POINT(goutin.coord),'bo');
patch('faces',coutout,'vertices',nodecoord,'facecolor','m','edgecolor','k');
plot(POINT(goutout.coord),'ko');

axis square
