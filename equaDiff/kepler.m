clear variables;
close all;

%paramètres
Grav = 4*pi^2;
T = 1;
pas = 0.0001;
tmin = 0;
tmax = 5*T;
%condition initiale
x0 = 0.5;
y0 = 0;
xp0 = 0;
yp0 = 11.5;

e=@(t,x,y,zx,zy) zx;
f=@(t,x,y,zx,zy) zy;
g=@(t,x,y,zx,zy) (-Grav*x)./((x.^2+y.^2).^(3/2));
h=@(t,x,y,zx,zy) (-Grav*y)./((x.^2+y.^2).^(3/2));

t= tmin:pas:tmax;

%Euler
xEuler = zeros(1,length(t));
yEuler = zeros(1,length(t));
zx = zeros(1,length(t));
zy =zeros(1,length(t));
xEuler(1) = x0;
yEuler(1) = y0;
zx(1) = xp0;
zy(1) = yp0;
for k=1:length(t)-1
  xEuler(k+1) = xEuler(k) + pas*e(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  yEuler(k+1) = yEuler(k) + pas*f(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  zx(k+1) = zx(k) + pas*g(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  zy(k+1) = zy(k) + pas*h(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
end

figure(1)
plot(xEuler, yEuler,0,0,'r*')

%Runge Kutta
xRK = zeros(1,length(t));
yRK = zeros(1,length(t));
zx = zeros(1,length(t));
zy =zeros(1,length(t));
xRK(1) = x0;
yRK(1) = y0;
zx(1) = xp0;
zy(1) = yp0;
for k=1:length(t)-1
  k1x = e(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  k1y = f(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  k1zx = g(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  k1zy = h(t(k), xEuler(k), yEuler(k), zx(k), zy(k));
  
  k2x = e(t(k) + pas/2, xEuler(k) + (pas/2)*k1x, yEuler(k) + (pas/2)*k1y, zx(k) + (pas/2)*k1zx, zy(k) + (pas/2)*k1zy);
  k2y = f(t(k) + pas/2, xEuler(k) + (pas/2)*k1x, yEuler(k) + (pas/2)*k1y, zx(k) + (pas/2)*k1zx, zy(k) + (pas/2)*k1zy);
  k2zx = g(t(k) + pas/2, xEuler(k) + (pas/2)*k1x, yEuler(k) + (pas/2)*k1y, zx(k) + (pas/2)*k1zx, zy(k) + (pas/2)*k1zy);
  k2zy = h(t(k) + pas/2, xEuler(k) + (pas/2)*k1x, yEuler(k) + (pas/2)*k1y, zx(k) + (pas/2)*k1zx, zy(k) + (pas/2)*k1zy);
  
  k3x = e(t(k) + pas/2, xEuler(k) + (pas/2)*k2x, yEuler(k) + (pas/2)*k2y, zx(k) + (pas/2)*k2zx, zy(k) + (pas/2)*k2zy);
  k3y = f(t(k) + pas/2, xEuler(k) + (pas/2)*k2x, yEuler(k) + (pas/2)*k2y, zx(k) + (pas/2)*k2zx, zy(k) + (pas/2)*k2zy);
  k3zx = g(t(k) + pas/2, xEuler(k) + (pas/2)*k2x, yEuler(k) + (pas/2)*k2y, zx(k) + (pas/2)*k2zx, zy(k) + (pas/2)*k2zy);
  k3zy = h(t(k) + pas/2, xEuler(k) + (pas/2)*k2x, yEuler(k) + (pas/2)*k2y, zx(k) + (pas/2)*k2zx, zy(k) + (pas/2)*k2zy);
  
  k4x = e(t(k) + pas, xEuler(k) + (pas)*k3x, yEuler(k) + (pas)*k3y, zx(k) + (pas)*k3zx, zy(k) + (pas)*k3zy);
  k4y = f(t(k) + pas, xEuler(k) + (pas)*k3x, yEuler(k) + (pas)*k3y, zx(k) + (pas)*k3zx, zy(k) + (pas)*k3zy);
  k4zx = g(t(k) + pas, xEuler(k) + (pas)*k3x, yEuler(k) + (pas)*k3y, zx(k) + (pas)*k3zx, zy(k) + (pas)*k3zy);
  k4zy = h(t(k) + pas, xEuler(k) + (pas)*k3x, yEuler(k) + (pas)*k3y, zx(k) + (pas)*k3zx, zy(k) + (pas)*k3zy);
  
  
  xRK(k+1) = xRK(k) + (pas/6)*(k1x+2*k2x+2*k3x+k4x);
  yRK(k+1) = yRK(k) + (pas/6)*(k1y+2*k2y+2*k3y+k4y);
  zx(k+1) = zx(k) + (pas/6)*(k1zx+2*k2zx+2*k3zx+k4zx);
  zy(k+1) = zy(k) + (pas/6)*(k1zy+2*k2zy+2*k3zy+k4zy);
  
end

figure(2)
plot(xRK, yRK,0,0,'r*')






