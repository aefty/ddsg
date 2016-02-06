close all;
xx = 0 : 0.01 : 1;

X = [1/8 1/4 ];
Y = [1 0];
[P,R,S] = lagrangepoly(X,Y);
plot(xx,polyval(P,xx),X,Y,'or',R,S,'.b',xx,spline(X,Y,xx),'--g'); hold on;


X = [1/4 3/8 1/2 ];
Y = [0 1 0];
[P,R,S] = lagrangepoly(X,Y);
plot(xx,polyval(P,xx),X,Y,'or',R,S,'.b',xx,spline(X,Y,xx),'--g'); hold on;

X = 1-[1/4 3/8 1/2 ];
Y = [0 1 0];
[P,R,S] = lagrangepoly(X,Y);
plot(xx,polyval(P,xx),X,Y,'or',R,S,'.b',xx,spline(X,Y,xx),'--g'); hold on;

axis([0 1 0 4])


grid
